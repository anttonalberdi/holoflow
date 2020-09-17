#17.09.2020 - Holoflow 0.1.

import subprocess
import argparse
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-i1', help="path1", dest="iread1", required=True)
parser.add_argument('-i2', help="path2", dest="iread2", required=True)
parser.add_argument('-bin_dir', help="input bin directory", dest="bin_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
#parser.add_argument('-R', help="Complete read group header line", dest="R", required=True)
args = parser.parse_args()

read1=args.read1
read2=args.read2
bin_dir=args.bin_dir
out_dir=args.out_dir
t=args.t
sample=args.sample
log=args.log
#R=args.R


# Run
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

    binlist = glob.glob(str(bin_dir)+"/*.fa")
    for bin in binlist:
        full_bin = os.path.abspath(bin)

            # define output files
        obam=''+out_dir+'/'+bin+'.bam'
        oread1=''+out_dir+'/'+bin+'_1.fastq'
        oread2=''+out_dir+'/'+bin+'_2.fastq'

    #Map bin to 1,2.fastq
        mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+full_bin+' '+iread1+' '+iread2+' | samtools view -T '+full_bin+' -b - > '+obam+''
        subprocess.check_call(mapCmd, shell=True)

        fastqCmd = 'module load tools samtools/1.9 && samtools view -T '+full_bin+' -b -f12 '+obam+' | samtools fastq -1 '+oread1+' -2 '+oread2+' -'
        subprocess.check_call(refbam2Cmd, shell=True)

        rmObamCmd = 'rm '+obam+''
        subprocess.check_call(rmAllbamCmd, shell=True)
