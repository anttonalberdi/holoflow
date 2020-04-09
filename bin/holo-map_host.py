#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-hostrg', help="host reference genome", dest="host_ref_gen", required=True)
parser.add_argument('-obam', help="all bam file", dest="all_bam", required=True)
args = parser.parse_args()

all_bam=args.all_bam
read1=args.read1
read2=args.read2
host_ref_gen=args.host_ref_gen

# Run
mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t 28 -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+host_ref_gen+' '+read1+' '+read2+' | samtools view -T '+host_ref_gen+' -b - > '+all_bam'
subprocess.check_call(mapCmd, shell=True)
