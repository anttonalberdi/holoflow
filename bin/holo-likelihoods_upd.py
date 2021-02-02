## 02.02.21 - Holoflow 0.1
import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bam_dir', help="bam files directory", dest="bam_dir", required=True)
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


    # Define extension of input files depending on variant caller
    if vc == "angsd":
        in_extension = '.beagle.gz'
    else:
        in_extension = '.vcf.gz'

    # Get all input files paths
    #in_variants = glob.glob(var_dir+'/*'+in_extension)


    # Run Beagle for chromosome
    for i in range(len(chr_list)):
        CHR = chr_list[i]
        #in_file = in_variants[i]
        in_file = 




module load java/1.8.0
module load bcftools/1.9
module load anaconda3/4.4.0

# update likelihoods
java -Xmxg -jar /services/tools/beagle/4.1/beagle.27Jul16.86a.jar gl=LD.vcf.gz ref=panel.vcf.gz chrom=${CHROM} gprobs=true out=${CHROM}_probs

bcftools index ${CHROM}_prob.vcf.gz

bcftools +setGT ${CHROM}_prob.vcf.gz -- -t q -n . -e'FORMAT/GP>=0.99' > ${CHROM}_prob_filt.vcf

bgzip ${CHROM}_prob_filt.vcf



#### - Input can be .vcf or beagle file : GATK,BCF//ANGSD
##### - gl= sample.chr name FILE
##### - ref=panel.vcf.gz is ref panel from filtering HD
##### bcftools -e'FORMAT/GP>=0.99' --> Those variants with likelihoods higher than 0.99, set as genotypes from which imputation the rest + refpanel
