
ANGSD:
module load htslib/1.9 angsd/0.931

angsd -bam sample_list.txt  -doGlf 2 -GL 1 -doPost 1 -doMaf 1 -doMajorMinor 1 -nThreads 10 -out file

parametros:
-GL con este parámetro se elige el modelo. 1 es para samtools. 2 para GATK. Estas dos opciones entiendo que son los que más nos interesan.
-doGLf outputs log genotype likehoods to a file.
-doMajorMinor 1 o 2. con 1 estima los major y minor alleles basandose en likelihoods data. Con la opción 2, a partir de recuentos de datos.
-doPost estimate posterior genotype probability based on the allele frequency as a prior
-doMaf frequency estimation. Opciones 1,2,4,8.
-nThreads



-out file name
    --> Snakefile specified

    
*no he adivinado todavía cómo definir el cromosoma.
http://www.popgen.dk/angsd/index.php/ANGSD
