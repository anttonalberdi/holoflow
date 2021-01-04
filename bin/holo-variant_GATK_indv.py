
GATK (es un poco más pesado):
module load java/1.8.0 gatk/4.1.8.1


- lista de BAM files. Por cada muestra una linea, tiene que aparecer todo el path de la muestra.
            --->  globglob

##################
Primero este paso hay que hacerlo para cada muestra por individual y por cromosoma: (GET SAMPLE ID FROM ARGPARSE)

for bam in bam_list: ##################

    bam_id = ...

    for chr in chr_list:

        gatk HaplotypeCaller --java-options "-XmxXXg" -R ${REF}  -I input.bam --ERC GVCF --native-pair-hmm-threads ${THREADS} --sample-ploidy 2 --min-prunning 1 --min-dangling-branch-length1 -L ${CHROM} -O ${BAM_ID}.raw.g.vcf.gz


        Estos parametros deberían ser opcionales, son para conseguir variantes más agresivos.
        --min-prunning 1
        --min-dangling-branch-length1
