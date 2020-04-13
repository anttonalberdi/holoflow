#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-i1', help="path1 input", dest="read1i", required=True)
parser.add_argument('-i2', help="path2 input", dest="read2i", required=True)
parser.add_argument('-o1', help="path1 output", dest="read1o", required=True)
parser.add_argument('-o2', help="path2 output", dest="read2o", required=True)
parser.add_argument('-a1', help="adapter 1 sequence", dest="a1", required=True)
parser.add_argument('-a2', help="adapter 2 sequence", dest="a2", required=True)
parser.add_argument('-maxns', help="max number of N's", dest="maxns", required=True)
parser.add_argument('-minq', help="minimum quality", dest="minq", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

read1i=args.read1i
read2i=args.read2i
read1o=args.read1o
read2o=args.read2o
a1=args.a1
a2=args.a2
maxns=args.maxns
minq=args.minq
threads=args.threads

# Run
qualfiltCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.1.3 && AdapterRemoval --file1 '+read1i+' --file2 '+read2i+' --output1 '+read1o+' --output2 '+read2o+' --trimqualities --trimns --maxns '+maxns+' --minquality '+minq+' --threads '+threads+' --adapter1 '+a1+' --adapter2'+a2+''
subprocess.check_call(qualfiltCmd, shell=True)
