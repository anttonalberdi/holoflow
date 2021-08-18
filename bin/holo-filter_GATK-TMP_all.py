## 26.02.21 - Holoflow 0.1
import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-var_dir', help="variant files directory", dest="var_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-QUAL', help="QUAL", dest="QUAL", required=True)
parser.add_argument('-QD', help="QD", dest="QD", required=True)
parser.add_argument('-FS', help="FS", dest="FS", required=True)
args = parser.parse_args()


var_dir=args.var_dir
out_dir=args.out_dir
chr_list=args.chr_list
ID=args.ID
log=args.log
threads=args.threads
QUAL=args.QUAL
QD=args.QD
FS=args.FS

## Run
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tFiltering of HD data with GATK - '+ID+'\n')
        logi.write(' \n\n')


    chromosome_list = list()
    # if the reference genome is not split by chromosomes but by scaffolds (for example)
    # remove -r region option and analyse all at once.
    # For this, chr_list will have only ONE row with 'ALL'
    all_genome_atonce = False
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            if chr.strip() == 'ALL':
                all_genome_atonce = True
            else:
                pass
            chromosome_list.append(chr.strip())


    for CHR in chromosome_list:
        geno_input = var_dir+'/'+ID+'.all_'+CHR+'.vcf'
        filter_output = out_dir+'/'+ID+'.HD_filt_'+CHR+'.vcf.gz'
        select_output = out_dir+'/'+ID+'.HD_filt_SNPs_'+CHR+'.vcf.gz'

        filterCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk VariantFiltration -V '+geno_input+' -filter "QD < '+QD+'" --filter-name "QD" -filter "QUAL < '+QUAL+'" --filter-name "QUAL" -filter "FS > '+FS+'" --filter-name "FS" -O '+filter_output+''
        subprocess.Popen(filterCmd,shell=True).wait()

        selectCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk SelectVariants -V '+filter_output+' --exclude-filtered --select-type-to-include SNP -O '+select_output+''
        subprocess.Popen(selectCmd,shell=True).wait()
