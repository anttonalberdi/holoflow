#############
BCFtools:
module load samtools/1.9 bcftools/1.9

samtools index ${SAMPLE}_map2host.bam
bcftools mpileup  -C 10 -q 10 -Q 10 -Ou -f ${REF} -r ${CHROM} -b sample_list.txt | bcftools call -m -v  -Oz -o all_${CHROM}.vcf.gz
bcftools view -m2 -M2 -v snps -Oz -o SNPs_${CHROM}.vcf.gz all_${CHROM}.vcf.gz

mpileup parameters:
-C coeficiente para degradar la calidad del mapeo. si se usa bwa, se recomienda usar 50
-q calidad de mapeo mínima
-Q calidad de base mínima
-r región, por cromosoma
-b lista de BAM files, en formato texto. Por cada muestra una linea, tiene que aparecer todo el path de la muestra.

call parameters:
-m multicaller mode
-v sólo llamar a variantes, no indels

view parameters: Este paso es para quedarse con los variantes bialélicos, sólo con snps.
http://samtools.github.io/bcftools/bcftools.html
