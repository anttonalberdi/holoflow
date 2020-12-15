
GATK (es un poco más pesado):
module load java/1.8.0 gatk/4.1.8.1

Primero este paso hay que hacerlo para cada muestra por individual y por cromosoma:
gatk HaplotypeCaller --java-options "-XmxXXg" -R ${REF}  -I input.bam --ERC GVCF --native-pair-hmm-threads ${THREADS} --sample-ploidy 2 --min-prunning 1 --min-dangling-branch-length1 -L ${CHROM} -O ${SAMPLE}.raw.g.vcf.gz


Estos parametros deberían ser opcionales, son para conseguir variantes más agresivos.
--min-prunning 1
--min-dangling-branch-length1

Después para todas las muestras a la vez por cromosoma:
gatk GenomicsDBImport --java-options "-Xmx XX g" --sample-name-map cohort.sample_map --genomicsdb-workspace-path my_database --reader-threads ${THREADS} -L ${CHROM}  -O ${SAMPLE}.raw.g.vcf.gz

gatk GenotypeGVCFs --java-options "-Xmx XX g" -R ${REF} -L ${CHROM}  -V gendb://my_database -O combined.raw.vcf

gatk GatherVcfs --java-options "-Xmx XX g" -I input -O output

gatk SelectVariants -V combined.raw.vcf  --select-type-to-include SNP -O SNPs_${CHROM} vcf.gz
#############
