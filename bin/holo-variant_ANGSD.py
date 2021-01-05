
ANGSD:
module load htslib/1.9 angsd/0.931


-b lista de BAM files, en formato lista? Por cada muestra una linea, tiene que aparecer todo el path de la muestra.
            1. --->  globglob
            2. write sample_list.txt file for file in globglob

-chr find out HOW TO SPECIFY CHR

-out_file = Snakefile_given_out_dir+group_name



            angsd -bam sample_list.txt  -doGlf 2 -GL 1 -doPost 1 -doMaf 1 -doMajorMinor 1 -nThreads 10 -out out_file











parametros:
-GL con este parámetro se elige el modelo. 1 es para samtools. 2 para GATK. Estas dos opciones entiendo que son los que más nos interesan.
-doGLf outputs log genotype likehoods to a file.
-doMajorMinor 1 o 2. con 1 estima los major y minor alleles basandose en likelihoods data. Con la opción 2, a partir de recuentos de datos.
-doPost estimate posterior genotype probability based on the allele frequency as a prior
-doMaf frequency estimation. Opciones 1,2,4,8.
-nThreads


######################################
######################################
######################################
IF LEARN HOW TO SPECIFY CHROMOSOME, LOOP OVER CHR LIST

-out file name
    --> Snakefile specified


*no he adivinado todavía cómo definir el cromosoma.
http://www.popgen.dk/angsd/index.php/ANGSD


######################################
