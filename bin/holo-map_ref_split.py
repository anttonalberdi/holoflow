#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-refg', help="reference genomes", dest="ref_gen", required=True)
parser.add_argument('-ibam', help="all bam file", dest="all_bam", required=True)
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-obam', help="bam file", dest="bam", required=True)
args = parser.parse_args()

all_bam=args.all_bam
ref_gen=args.ref_gen
bam=args.bam
read1=args.read1
read2=args.read2

# Run
refbam1Cmd = 'module load tools samtools/1.9 && samtools view -T '+ref_gen+' -b -F12 '+all_bam+' > '+bam+''
subprocess.check_call(refbam1Cmd, shell=True)

refbam2Cmd = 'module load tools samtools/1.9 && samtools view -T '+ref_gen+' -b -f12 '+all_bam+' | samtools fastq -1 '+read1+' -2 '+read2+' -'
subprocess.check_call(refbam2Cmd, shell=True)

rmAllbamCmd = 'rm '+all_bam+''
subprocess.check_call(rmAllbamCmd, shell=True)
