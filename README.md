# Quick phylogeny using BUSCO

#### Nextflow pipeline to use BUSCO results (from DNA) to create a phylogenetic tree

## 1 Tools and dependencies
Tools used from other projects are shown below.

- BUSCO
- Augustus
- Seqtk
- Mafft
- Gblocks
- FASconCAT
- ProtTest

Then, clone the repository:

`git clone https://github.com/AllanSSX/quick_phylo_BUSCO.git`

And run the pipelines

`nextflow run phylo-busco.nf *.fna`

The results are located in the `results` folder. The 2 importants output files are:

- `FcC_smatrix.fas` : fasta file containing the alignments. Input file for RAxML
- `prottest.txt` : Prottest result. Allows you to choose the correct model for RAxML

