# holoflow
Bioinformatics pipeline for hologenomics data generation and analysis

module unload gcc/5.1.0
module load anaconda3/4.4.0

snakemake -s Snakefile -n -r ${workdir}/02-DuplicatesRemoved/H2A_1.fastq ${workdir}/02-DuplicatesRemoved/H2A_2.fastq
