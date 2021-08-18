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
parser.add_argument('-QUAL', help="QUAL", dest="QUAL", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


var_dir=args.var_dir
out_dir=args.out_dir
chr_list=args.chr_list
QUAL=args.QUAL
ID=args.ID
log=args.log
threads=args.threads


## Run
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tFiltering of HD data with BCFtools - '+ID+'\n')
        logi.write(' \n\n')

    chromosome_list = list()
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            chromosome_list.append(chr.strip())

    for CHR in chromosome_list:
        mpileup_input = var_dir+'/'+ID+'.all_'+CHR+'.vcf.gz'
        filter_output = out_dir+'/'+ID+'.HD_filt_'+CHR+'.vcf.gz'
        view_output = out_dir+'/'+ID+'.HD_SNPs_'+CHR+'.vcf.gz'

        # Filter variants by quality and depth
        filterCmd='module load bcftools/1.11 && bcftools filter -s LowQual -e "%QUAL<'+QUAL+' || DP<(AVG(DP)*3)" --threads '+threads+' -Oz -o '+filter_output+' '+mpileup_input+''
        subprocess.Popen(filterCmd,shell=True).wait()

        viewCmd='module load bcftools/1.11 && bcftools view -m2 -M2 -v snps --threads '+threads+' -Oz -o '+view_output+' '+filter_output+''
        subprocess.Popen(viewCmd,shell=True).wait()
