
params.fasta = "*.fna"
params.outdir = "results"

genomesFasta = Channel
	.fromPath( params.fasta )
	.map { file -> tuple(file.baseName, file) }

/*
 * RUN ANALYSIS
 */
 
process busco {
	
	module 'busco/3.1.0'
	
	publishDir "${params.outdir}/1-busco", mode: 'copy'
	
	input:
	set name, file(fasta) from genomesFasta
	
	output:
	file ("run_${name}/*.{txt,tsv}") into busco_summary_results
	file ("run_${name}/single_copy_busco_sequences_${name}.faa") into busco_single_copy_proteins
	
	script:
	
	"""
	export AUGUSTUS_CONFIG_PATH=/home/acormier/soft/Augustus/3.3.2/config
	run_BUSCO.py -i ${fasta} -o ${name} -m genome -l /home/databases/busco/microsporidia_odb9/ --cpu 6
	cat run_${name}/single_copy_busco_sequences/*.faa > run_${name}/single_copy_busco_sequences_${name}.faa
	"""
}

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

process filter_single_copy {
	
	publishDir "${params.outdir}/2-single-copy", mode: 'copy'
	
	input:
	file fasta from busco_single_copy_proteins_concat
	
	output:
	file "*.lst" into busco_single_copy_proteins_list mode flatten
	
	// add a limitation on -s option in order than it can be superior to the number of input fasta
	script:
	"""
	single_copies_busco.py -f ${fasta} -s 3
	"""		
}

process seqtk {
	
	module 'seqtk/1.3-r106'
	
	publishDir "${params.outdir}/2-single-copy", mode: 'copy'
	
	input:
	//file lst from busco_single_copy_proteins_list.collect()
	file lst from busco_single_copy_proteins_list
	file fasta from busco_single_copy_proteins_concat

	output:
	file "*.faa" into busco_single_copy_proteins_toMafft

	script:
	"""
	seqtk subseq -l 70 ${fasta} ${lst} > ${lst}.faa
	"""
	//for liste in ${lst}
	//do
	//	seqtk subseq -l 70 ${fasta} \${liste} > \${liste}.faa
	//done
}

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
	"""
}

process gblocks {
	
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

process cleanHeader {
	
	module 'Gblocks/0.91b'
	
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

process FASconCAT {
	
	module 'FASconCAT/1.0'
	
	publishDir "${params.outdir}/5-cat", mode: 'copy'
	
	input:
	file gblocks from busco_single_copy_proteins_toFASconCAT
	
	output:
	file "FcC_smatrix.fas" into busco_single_copy_proteins_toProTest
	
	script:
	"""
	FASconCAT_v1.0.pl -s
	"""
}


process ProtTest {
	
	module 'ProtTest/3.4.2'
	
	publishDir "${params.outdir}/5-cat", mode: 'copy'
	
	input:
	file matrix from busco_single_copy_proteins_toProTest
	
	output:
	file "prottest.txt" into busco_single_copy_proteins_toProTestResult
	
	script:
	"""
	prottest3 -i ${matrix} -all-distributions -F -AIC -BIC -tc 0.5 -o prottest.txt -threads 4
	"""
}


