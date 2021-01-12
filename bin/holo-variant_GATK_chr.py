## 11.01.20 - Holoflow 0.1

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bam_dir', help="bam files directory", dest="bam_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ref_g', help="reference genome", dest="ref_g", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-degr_mapp_qual', help="degradation mapping quality", dest="degr_mqual", required=True)
parser.add_argument('-min_mapp_qual', help="minimum mapping quality", dest="min_mqual", required=True)
parser.add_argument('-min_base_qual', help="minimum base quality", dest="min_bqual", required=True)
parser.add_argument('-chr_region', help="specific chromosome region", dest="chr_region", required=True)
parser.add_argument('-multicaller', help="multicaller option", dest="multicaller", required=True)
#parser.add_argument('-not_indels', help="only variants not indels", dest="not_indels", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


bam_dir=args.bam_dir
out_dir=args.out_dir
ref_g=args.ref_g
chr_list=args.chr_list
degr_mqual=args.degr_mqual
min_mqual=args.min_mqual
min_bqual=args.min_bqual
chr_region=args.chr_region
multicaller=args.multicaller
#not_indels=args.not_indels
ID=args.ID
log=args.log
threads=args.threads

## Run
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tVariant calling with BCFtools tep - '+ID+'\n')
        logi.write(' \n\n')

    # Get chromosomes list
    chromosome_list = list()
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            chromosome_list.append(chr.strip())



##############################
GATK (es un poco más pesado):
module load java/1.8.0 gatk/4.1.8.1


- lista de BAM files. Por cada muestra una linea, tiene que aparecer todo el path de la muestra.
            --->  globglob


##################
Después para todas las muestras a la vez por cromosoma: (ID)
for chr in chr_list: ##################


### Isn't GenomicsDBImport supposed to go before this chr loop? inside the by-sample loop
    gatk GenomicsDBImport --java-options "-Xmx28g" --sample-name-map cohort.sample_map --genomicsdb-workspace-path ${PATH OUT}/my_database --reader-threads ${THREADS} -L ${CHR} 2> >(tee "$logfile")

    gatk GenotypeGVCFs --java-options "-Xmx XX g" -R ${REF} -L ${CHROM}  -V gendb://my_database -O combined.raw.vcf

    gatk GatherVcfs --java-options "-Xmx XX g" -I input -O output

    gatk SelectVariants -V combined.raw.vcf  --select-type-to-include SNP -O SNPs_${CHROM}.vcf.gz
    #############
