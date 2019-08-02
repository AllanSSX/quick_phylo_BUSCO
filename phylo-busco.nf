//------------------------------------------
// PIPELINE INPUT PARAMETERS 
//------------------------------------------

import ParamsHelper

params.indir = "."
params.outdir = "results"
ParamsHelper.checkNonEmptyParam(params.species, "species");

BUSCOodb = params.odb
augustusConfig = params.augustus

//------------------------------------------
// RUN ANALYSIS
//------------------------------------------

// Create a channel for all fna (genomes) files
genomesFasta = Channel
	.fromPath( [params.indir, '*.fna'].join(File.separator) )
	.map { file -> tuple(file.baseName, file) }


// Launch BUSCO and cat all single copy genes into a single file for each specie
process busco {
	cpus params.busco.cpus
	module 'busco/3.1.0'
	
	publishDir "${params.outdir}/1-busco", mode: 'copy'
	
	input:
	set name, file(fasta) from genomesFasta
	augustusConfig
	BUSCOodb
	
	output:
	file ("run_${name}/*.{txt,tsv}") into busco_summary_results
	file ("run_${name}/single_copy_busco_sequences_${name}.faa") into busco_single_copy_proteins
	
	script:
	
	"""
	export AUGUSTUS_CONFIG_PATH=${augustusConfig}
	run_BUSCO.py -i ${fasta} -o ${name} -m genome -l ${BUSCOodb} --cpu ${params.busco.cpus}
	cat run_${name}/single_copy_busco_sequences/*.faa > run_${name}/single_copy_busco_sequences_${name}.faa
	"""
}

// Concat all single copy genes of all specie into a single file
process concat_busco {
	
	publishDir "${params.outdir}/1-busco", mode: 'copy'
	
	input:
	file fasta from busco_single_copy_proteins.collect()
	
	output:
	file "single_copy_busco_sequences.fasta" into busco_single_copy_proteins_concat
	
	script:
	"""
	cat ${fasta} > single_copy_busco_sequences.fasta
	"""
}

// Filter in order to keep at least a minimum of X species which shared a single copy gene (give a list of sequence ID)
process filter_single_copy {
	
	//module 'anaconda3/conda3'
	
	publishDir "${params.outdir}/2-single-copy", mode: 'copy'
	
	input:
	file fasta from busco_single_copy_proteins_concat
	
	output:
	file "*.lst" into busco_single_copy_proteins_list mode flatten
	// flatten to allow parallelization at the next step
	
	// add a limitation on -s option in order than it can be superior to the number of input fasta
	script:
	"""
	single_copies_busco.py -f ${fasta} -s params.species
	"""		
}

// Extract each ortholgous sequences
process seqtk {
	
	module 'seqtk/1.3-r106'
	
	publishDir "${params.outdir}/2-single-copy", mode: 'copy'
	
	input:
	file lst from busco_single_copy_proteins_list
	file fasta from busco_single_copy_proteins_concat

	output:
	file "*.faa" into busco_single_copy_proteins_toMafft

	script:
	"""
	seqtk subseq -l 70 ${fasta} ${lst} > ${lst}.faa
	"""
}

// Align each orthogroup
process mafft {
	
	module 'mafft/7.407'
	
	publishDir "${params.outdir}/3-mafft", mode: 'copy'
	
	input:
	file faa from busco_single_copy_proteins_toMafft
	
	output:
	file "*.mafft" into busco_single_copy_proteins_toGblocks
	
	script:
	"""
	mafft --auto ${faa} > ${faa}.mafft
	sleep 4
	"""
}

// Clean the alignment
process gblocks {
	// As Gblocks exit status is always 1...
	validExitStatus 1
	
	module 'Gblocks/0.91b'
	
	publishDir "${params.outdir}/4-gblocks", mode: 'copy'
	
	input:
	file aln from busco_single_copy_proteins_toGblocks
	
	output:
	file "*-gb" into busco_single_copy_proteins_toCleanHeader
	
	script:
	"""
	Gblocks ${aln} -t=p -p=n -b3=8 -b4=10 -b5=h
	"""
}

// Clean the header for higher readability in the final tree
// Default: EOG093001PG:NGRA.fna:NGRA000074:9808-10725
// Cleaned: NGRA
process cleanHeader {
	
	publishDir "${params.outdir}/5-cat", mode: 'copy'
	
	input:
	file gblocks from busco_single_copy_proteins_toCleanHeader
	
	output:
	file "*.fas" into busco_single_copy_proteins_toFASconCAT
	
	script:
	"""
	cleanHeaderFromBuscoSG.py -f ${gblocks} > ${gblocks}.fas
	"""
}

// Create the matrix for RAxML and ProtTest
process FASconCAT {
	
	module 'FASconCAT/1.0'
	
	publishDir "${params.outdir}/5-cat", mode: 'copy'
	
	input:
	file gblocks from busco_single_copy_proteins_toFASconCAT.collect()
	
	output:
	file "FcC_smatrix.fas" into busco_single_copy_proteins_toProTest
	
	script:
	"""
	FASconCAT_v1.0.pl -s
	"""
}

// Find the better paramters for RAxML
process ProtTest {
	cpus params.prottest.cpus
	module 'ProtTest/3.4.2'
	
	publishDir "${params.outdir}/5-cat", mode: 'copy'
	
	input:
	file matrix from busco_single_copy_proteins_toProTest
	
	output:
	file "prottest.txt" into busco_single_copy_proteins_toProTestResult
	
	script:
	"""
	prottest3 -i ${matrix} -all-distributions -F -AIC -BIC -tc 0.5 -o prottest.txt -threads ${params.prottest.cpus}
	"""
}