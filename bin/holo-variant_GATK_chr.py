
GATK (es un poco más pesado):
module load java/1.8.0 gatk/4.1.8.1


- lista de BAM files. Por cada muestra una linea, tiene que aparecer todo el path de la muestra.
            --->  globglob


##################
Después para todas las muestras a la vez por cromosoma: (ID)
for chr in chr_list: ##################


### Isn't GenomicsDBImport supposed to go before this chr loop? inside the by-sample loop
    gatk GenomicsDBImport --java-options "-Xmx XX g" --sample-name-map cohort.sample_map --genomicsdb-workspace-path my_database --reader-threads ${THREADS} -L ${CHROM}  -O ${SAMPLE}.raw.g.vcf.gz

    gatk GenotypeGVCFs --java-options "-Xmx XX g" -R ${REF} -L ${CHROM}  -V gendb://my_database -O combined.raw.vcf

    gatk GatherVcfs --java-options "-Xmx XX g" -I input -O output

    gatk SelectVariants -V combined.raw.vcf  --select-type-to-include SNP -O SNPs_${CHROM}.vcf.gz
    #############
