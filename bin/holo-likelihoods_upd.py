## 02.02.21 - Holoflow 0.1
import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-var_dir', help="variant files directory", dest="var_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ref_panel', help="reference panel", dest="ref_panel", required=True)
parser.add_argument('-vc', help="variant caller", dest="vc", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


var_dir=args.var_dir
out_dir=args.out_dir
ref_panel=args.ref_panel
vc=args.vc
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
        logi.write('\t\t'+current_time+'\tLikelihoods update with Beagle for Low Depth samples step - '+ID+'\n')
        logi.write(' \n\n')

    # Get file extension depending on variant caller
    if vc == "angsd":
        in_extension = '.beagle.gz'
    else:
        in_extension = '.vcf.gz'


    # Run Beagle per chromosome
    for CHR in chr_list:

        in_file_base = var_dir+'/'+ID+'.SNPs_'+CHR+in_extension
        bgl_out_base = out_dir+'/'+ID+'.probs_'+CHR

        bglCmd = 'module load java/1.8.0 anaconda3/4.4.0 && java  -jar /services/tools/beagle/4.1/beagle.27Jul16.86a.jar gl='+in_file_base+' ref='+ref_panel+' chrom='+CHR+' gprobs=true out='+bgl_out_base+' t='+threads+''
        subprocess.Popen(bglCmd,shell=True).wait()

        # Index and set genotypes in output
        bgl_out = bgl_out_base+'.vcf.gz'
        filt_out = out_dir+'/'+ID+'.probs_filt.vcf'

        bcfCmd = 'module load tools bcftools/1.11 && bcftools index '+bgl_out+' && bcftools +setGT '+bgl_out+' -- -t -q -n . -e "FORMAT/GP>=0.99" > '+filt_out+' && bgzip '+filt_out+''
        subprocess.Popen(bcfCmd,shell=True).wait()
