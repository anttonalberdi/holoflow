## 11.01.20 - Holoflow 0.1

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-vcf_dir', help="individual vcf files directory", dest="vcf_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ref_g', help="reference genome", dest="ref_g", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-Dquality', help="data quality", dest="Dquality", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


vcf_dir=args.vcf_dir
out_dir=args.out_dir
ref_g=args.ref_g
chr_list=args.chr_list
Dquality=args.Dquality
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

    # Run GATK
    for CHR in chromosome_list:
        sample_map_name = vcf_dir+'/sample_map.'+CHR+'.txt'

        # Define outputs
        my_database = out_dir+'/'+CHR+'_database'
        geno_output = out_dir+'/'+ID+'.combined_'+CHR+'.raw.vcf'
        variants_output = out_dir+'/'+ID+'.LD_SNPs_'+CHR+'.vcf.gz'


        dbCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk GenomicsDBImport --java-options "-Xmx180g" --sample-name-map '+sample_map_name+' --genomicsdb-workspace-path '+my_database+' --reader-threads '+threads+' -L '+CHR+''
        subprocess.Popen(dbCmd,shell=True).wait()

        # -V gendb://my_database
        genoCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk GenotypeGVCFs --java-options "-Xmx180g" -R '+ref_g+' -L '+CHR+' -V gendb://'+my_database+' -O '+geno_output+''
        subprocess.Popen(genoCmd,shell=True).wait()

        if Dquality == 'LD':
            variantsCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk SelectVariants -V '+geno_output+'  --select-type-to-include SNP -O '+variants_output+''
            subprocess.Popen(variantsCmd,shell=True).wait()
        else:
            pass

        if CHR == chromosome_list[-1]:
            rmCmd='rm -rf '+vcf_dir+''
            subprocess.Popen(rmCmd,shell=True).wait()
