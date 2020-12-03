 #13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import re
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-fq_path', help="path to .fastq files", dest="fq_path", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-obam_b', help="output bam file base", dest="obam_base", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
fq_path=args.fq_path
t=args.t
obam_base=args.obam_base
ID=args.ID
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tCoAssembly Mapping step - '+ID+'\n')
    log.write('The original metagenomic reads are being mapped to the indexed assembly so coverage info can be retrieved.\n\n')


# Get read1 and read2 paths

reads1=glob.glob(fq_path+'/*_1.f*')

for read1 in reads:
    read1=os.path.basename(read1)
    sampleID=re.sub('(\.|_)[0-9]{1}\.f[aA-zZ]*\.?.*','',read1) # remove .1.fa .1.fastq _1.fq.gz _1.fastq.gz ...

    obam=obam_b+'/'+sampleID+'.mapped.bam'
    read2= read1.replace('1','2')

    if not os.path.exists(str(obam)):
        mappingCmd='module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+a+' '+fq_path+'/'+read1+' '+fq_path+'/'+read2+' | samtools view -b - | samtools sort -T '+sampleID+' -o '+obam+''
        subprocess.check_call(mappingCmd, shell=True)
