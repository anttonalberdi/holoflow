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
        logi.write('\t\t'+current_time+'\tVariant calling with BCFtools step - '+ID+'\n')
        logi.write(' \n\n')

    # Get chromosomes list
    chromosome_list = list()
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            chromosome_list.append(chr.strip())



    # Generate bam files' paths file list & index
    bam_list = [os.path.basename(x) for x in glob.glob(bam_dir+'/*.bam')]
    bam_list_file = out_path+'/'+ID+'_bam_list.txt'

    with open(bam_list_file,'w+') as bam_files:

        for bam in bam_list:
            bam_files.write(str(bam)+'\n')

            if not os.path.isfile(bam+'.bai'): # If not indexed, index bam - Theoretically these are sorted from preprocessing
                idxbamCmd = 'module load tools samtools/1.9 && samtools index '+bam+''
                subprocess.Popen(idxbamCmd,shell=True).wait()

            else:
                pass

    # Run BCFtools
    for CHR in chromosome_list:

        mpileup_output = out_dir+'/'+ID+'.all_'+CHR+'.vcf.gz'
        view_output = out_dir+'/'+ID+'.SNPs_'+CHR+'.vcf.gz'

        if not (chr_region == 'False'):

            if not (multicaller == 'False'):
                bcf1Cmd = 'module load tools samtools/1.9 bcftools/1.9 && bcftools mileup -C '+degr_mqual+' -q '+min_mqual+' -Q '+min_bqual+' -Ou  -f '+ref_g+' -r '+CHR+' -b '+bam_list_file+' -r '+chr_region+' | bcftools call -m -v -Oz -o '+mpileup_output+''
                subprocess.Popen(bcf1Cmd,shell=True).wait()
                bcf2Cmd = 'bcftools view -m2 -M2 -v snps -Oz -o '+view_output+' '+mpileup_output+''
                subprocess.Popen(bcf2Cmd,shell=True).wait()

            else:
                bcf1Cmd = 'module load tools samtools/1.9 bcftools/1.9 && bcftools mileup -C '+degr_mqual+' -q '+min_mqual+' -Q '+min_bqual+' -Ou  -f '+ref_g+' -r '+CHR+' -b '+bam_list_file+' -r '+chr_region+' | bcftools call -v -Oz -o '+mpileup_output+''
                subprocess.Popen(bcf1Cmd,shell=True).wait()
                bcf2Cmd = 'bcftools view -m2 -M2 -v snps -Oz -o '+view_output+' '+mpileup_output+''
                subprocess.Popen(bcf2Cmd,shell=True).wait()


        else:
            if not (multicaller == 'False'):
                bcf1Cmd = 'module load tools samtools/1.9 bcftools/1.9 && bcftools mileup -C '+degr_mqual+' -q '+min_mqual+' -Q '+min_bqual+' -Ou  -f '+ref_g+' -r '+CHR+' -b '+bam_list_file+' | bcftools call -m -v -Oz -o '+mpileup_output+''
                subprocess.Popen(bcf1Cmd,shell=True).wait()
                bcf2Cmd = 'bcftools view -m2 -M2 -v snps -Oz -o '+view_output+' '+mpileup_output+''
                subprocess.Popen(bcf2Cmd,shell=True).wait()

            else:
                bcf1Cmd = 'module load tools samtools/1.9 bcftools/1.9 && bcftools mileup -C '+degr_mqual+' -q '+min_mqual+' -Q '+min_bqual+' -Ou  -f '+ref_g+' -r '+CHR+' -b '+bam_list_file+' | bcftools call -v -Oz -o '+mpileup_output+''
                subprocess.Popen(bcf1Cmd,shell=True).wait()
                bcf2Cmd = 'bcftools view -m2 -M2 -v snps -Oz -o '+view_output+' '+mpileup_output+''
                subprocess.Popen(bcf2Cmd,shell=True).wait()
