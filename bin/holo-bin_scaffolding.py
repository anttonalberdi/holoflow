#10.09.2020 - Holoflow 0.1.

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



dt_bd=args.dt_bd
out_dir=args.out_dir
sample=args.sample
log=args.log
threads=args.threads


# Run
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

    #Create library file
        # Insertion size between paired reads: 150
        # Maximum allowed error: 1
    libCmd='mkdir '+out_dir+' && printf "'+sample+' bwa '+fq_dir+'/'+sample+'_1.fastq '+fq_dir+'/'+sample+'_2.fastq 150 1 FR" > '+out_dir+'/'+sample+'.lib'
    subprocess.check_call(libCmd, shell=True)


    #Run SSPACE
    binlist = glob.glob(str(bin_dir)+"/*.fa")
    for bin in binlist:
        full_bin = os.path.abspath(bin)
        bin_id = bin.replace(".contigs.fa","")
        sspaceCmd = 'cd '+outdir+' && cd .. && module load tools perl/5.24.0 sspace-standard/3.0 parallel/20190522 && SSPACE_Standard_v3.0.pl -l '+out_dir+'/'+sample+'.lib -s '+full_bin+' -x 1 -T '+threads+' -o 5 -m 16 -k 2 -n 10 -b '+bin_id+''
        subprocess.check_call(sspaceCmd, shell=True)


    #Rearrange outputs
    for bin in binlist:
        bin_id = bin.replace(".contigs.fa","")
        faoutpCmd='cp 'out_dir'/'+bin_id+'.final.scaffolds.fasta 'out_dir'/../'+bin_id+'.fa'
        subprocess.check_call(faoutpCmd, shell=True)
        infoutCmd='cp 'out_dir'/'+bin_id+'.summaryfile.txt 'out_dir'/../'+bin_id+'.info'
        subprocess.check_call(infoutCmd, shell=True)
    ##  rmCmd='rm 'out_dir''
    ##  subprocess.check_call(rmCmd, shell=True)
