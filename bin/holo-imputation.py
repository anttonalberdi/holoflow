## 02.02.21 - Holoflow 0.1
import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-upd_dir', help="updated likelihoods files directory", dest="upd_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ref_panel', help="reference panel", dest="ref_panel", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


upd_dir=args.upd_dir
out_dir=args.out_dir
ref_panel=args.ref_panel
chr_list=args.chr_list
ID=args.ID
log=args.log
threads=args.threads


## Run
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tGenotypes are being imputed using updated likelihoods with Beagle for Low Depth samples step - '+ID+'\n')
        logi.write(' \n\n')

    chromosome_list = list()
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            chromosome_list.append(chr.strip())

    for CHR in chromosome_list:

        in_file = upd_dir+'/'+ID+'.probs_'+CHR+'.vcf.gz'
        bgl_out_base = out_dir+'/'+ID+'.imputed_'+CHR

        # Run imputation

        bglCmd = 'module load java/1.8.0 anaconda3/4.4.0 && java -Xmx180g -jar /services/tools/beagle/5.1/beagle-5.1.jar gt='+in_file+' ref='+ref_panel+' chrom='+CHR+' gp=true out='+bgl_out_base+''
        subprocess.Popen(bglCmd,shell=True).wait()

        bgl_out = bgl_out_base+'.vcf.gz'
        bcf_out = out_dir+'/'+ID+'.imputed_filt_'+CHR+'.vcf'

        bcfCmd = 'module load bcftools/1.11 && bcftools index '+bgl_out+' && bcftools +setGT '+bgl_out+' -- -t q -n . -e"FORMAT/GP>=0.99" > '+bcf_out+' && bgzip '+bcf_out+''
        subprocess.Popen(bcfCmd,shell=True).wait()
