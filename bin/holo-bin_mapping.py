#17.09.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-i1', help="path1", dest="read1", required=True)
parser.add_argument('-i2', help="path2", dest="read2", required=True)
parser.add_argument('-bin_dir', help="input bin directory", dest="bin_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
#parser.add_argument('-R', help="Complete read group header line", dest="R", required=True)
args = parser.parse_args()

read1=args.read1
read2=args.read2
bin_dir=args.bin_dir
out_dir=args.out_dir
t=args.t
ID=args.ID
log=args.log
#R=args.R


# Run
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tBin Mapping step - '+ID+'\n')
        logi.write('This step retrieves the paired-end reads found in each bin as they are to be used in the next step.\n\n')


    binlist = glob.glob(str(bin_dir)+"/dereplicated_genomes/*.fa")
    for bin in binlist:
        bin_name=os.path.basename(bin)
        bin_name=bin_name.replace(".fa","")


            # define output files
        obam=''+out_dir+'/'+bin_name+'.bam'
        oread1=''+out_dir+'/'+bin_name+'_1.fastq'
        oread2=''+out_dir+'/'+bin_name+'_2.fastq'

    #Map bin to 1,2.fastq

        idxbwaCmd='module load tools bwa/0.7.15 && bwa index '+bin+''
        subprocess.check_call(idxbwaCmd, shell=True)

        idxsamCmd='module load tools samtools/1.11 && samtools faidx '+bin+''
        subprocess.check_call(idxsamCmd, shell=True)


        mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+bin+' '+read1+' '+read2+' | samtools view -T '+bin+' -b - > '+obam+''
        subprocess.check_call(mapCmd, shell=True)

        fastqCmd = 'module load tools samtools/1.11 && samtools view -T '+bin+' -b -f12 '+obam+' | samtools fastq -1 '+oread1+' -2 '+oread2+' -'
        subprocess.check_call(fastqCmd, shell=True)

        rmvbamCmd = 'rm '+obam+' '+bin+'.*'
        subprocess.check_call(rmvbamCmd, shell=True)
