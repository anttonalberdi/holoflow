#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-hostrg', help="host reference genome", dest="host_ref_gen", required=True)
parser.add_argument('-ibam', help="all bam file", dest="all_bam", required=True)
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-obam', help="host bam file", dest="host_bam", required=True)
args = parser.parse_args()

all_bam=args.all_bam
host_ref_gen=args.host_ref_gen
host_bam=args.host_bam
read1=args.read1
read2=args.read2

# Run
hostbam1Cmd = 'module load tools samtools/1.9 && samtools view -T '+host_ref_gen+' -b -F12 '+all_bam+' > '+host_bam+''
subprocess.check_call(hostbam1Cmd, shell=True)
hostbam2Cmd = 'module load tools samtools/1.9 && samtools view -T '+host_ref_gen+' -b -f12 '+all_bam+' | samtools fastq -1 '+read1+' -2 '+read2+' -'
subprocess.check_call(hostbam2Cmd, shell=True)
rmAllbamCmd = 'rm '+all_bam+''
subprocess.check_call(rmAllbamCmd, shell=True)
