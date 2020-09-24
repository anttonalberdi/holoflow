#17.09.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-fq_dir', help="input .fq directory", dest="fq_dir", required=True)
parser.add_argument('-bin_dir', help="input bin directory", dest="bin_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


fq_dir=args.fq_dir
bin_dir=args.bin_dir
out_dir=args.out_dir
sample=args.sample
log=args.log
threads=args.threads


# Run
# if not (os.path.exists(str(out_dir))):
#     os.mkdir(str(out_dir))

binlist = glob.glob(str(bin_dir)+"/*.fa")
for bin in binlist:
    bin_name=os.path.basename(bin)
    bin_name = bin_name.replace(".fa","")
    print(bin)
    print(bin_name)
    lib_file=str(out_dir+'/'+bin_name+'.lib')

    #Create library file
        # Insertion size between paired reads: 150
        # Maximum allowed error: 1
    libCmd='printf "'+sample+' bwa '+fq_dir+'/'+bin_name+'_1.fastq '+fq_dir+'/'+bin_name+'_2.fastq 150 1 FR" >> '+lib_file+''
    subprocess.check_call(libCmd, shell=True)

    #Run SSPACE
    sspaceCmd ='cd '+out_dir+' && module load tools perl/5.24.0 sspace-standard/3.0 parallel/20190522 && SSPACE_Standard_v3.0.pl -l '+lib_file+' -s '+bin+' -x 1 -T '+threads+' -o 5 -m 16 -k 2 -n 10 -b '+bin_name+''
    subprocess.check_call(sspaceCmd, shell=True)


#Rearrange outputs
# for bin in binlist:
#     bin_name = bin.replace(".contigs.fa","")
#     faoutpCmd='cp 'out_dir'/'+bin_name+'.final.scaffolds.fasta 'out_dir'/../'+bin_name+'.fa'
#     subprocess.check_call(faoutpCmd, shell=True)
#     infoutCmd='cp 'out_dir'/'+bin_name+'.summaryfile.txt 'out_dir'/../'+bin_name+'.info'
#     subprocess.check_call(infoutCmd, shell=True)
##  rmCmd='rm 'out_dir''
##  subprocess.check_call(rmCmd, shell=True)
