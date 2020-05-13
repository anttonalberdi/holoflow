#13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-1', help="read1", dest="read1", required=True)
parser.add_argument('-2', help="read2", dest="read2", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-obam', help="output bam file", dest="obam", required=True)
args = parser.parse_args()


a=args.a
read1=args.read1
read2=args.read2
t=args.t
obam=args.obam



# Run
mappingCmd='module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+a+' '+read1+' '+read2+' | samtools view -T '+a+' -b - | samtools sort -T '+a+' - > '+obam+''
subprocess.check_call(mappingCmd, shell=True)
