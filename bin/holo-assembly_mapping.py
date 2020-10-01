 #13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-1', help="read1", dest="read1", required=True)
parser.add_argument('-2', help="read2", dest="read2", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-obam', help="output bam file", dest="obam", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
read1=args.read1
read2=args.read2
t=args.t
obam=args.obam
sample=args.sample
log=args.log



# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tAssembly Mapping step - Sample '+sample+'\n')
    log.write('The original metagenomic reads are being mapped to the indexed assembly so coverage info can be retrieved.\n\n')


if not os.path.exists(str(obam)):
    mappingCmd='module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+a+' '+read1+' '+read2+' | samtools view -T '+a+' -b - | samtools sort -T '+a+' - > '+obam+''
    subprocess.check_call(mappingCmd, shell=True)
