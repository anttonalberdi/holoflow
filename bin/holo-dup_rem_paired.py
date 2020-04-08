#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-o', help="output directory", dest="output_dir", required=True)
parser.add_argument('-sep', help="sep", dest="separator", required=True)
args = parser.parse_args()

outpur_dir=args.output_dir
read1=args.read1
read2=args.read2
separator=args.separator

# Run
seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d'+separator+' '+read1+' '+read2+' | seqkit rmdup -s -j 28 -o'+ output_dir'
subprocess.check_call(seqkitCmd, shell=True)
