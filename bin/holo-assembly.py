#09.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-o', help="output directory", dest="out", required=True)
parser.add_argument('-m', help="memory", dest="memory", required=True)
parser.add_argument('-k_megahit', help="k-mer size list megahit", dest="k_megahit", required=True)
parser.add_argument('-k_spades', help="k-mer size list spades", dest="k_spades", required=True)
parser.add_argument('-a', help="assembler", dest="assembler", required=True)
args = parser.parse_args()


read1=args.read1
read2=args.read2
out=args.out
memory=args.memory
k_megahit=args.k_megahit
k_spades=args.k_spades
threads=args.threads
assembler=args.assembler

# Run
if assembler == "megahit":
    megahitCmd = shell('module load tools megahit/1.1.1 && megahit -1 '+read1+' -2 '+read2+' -t '+threads+' --k-list '+k_megahit+' -o '+out+'')
    subprocess.check_call(megahitCmd, shell=True)

if assembler == "spades":
    spadesCmd = shell('module unload anaconda3/4.4.0 && module load tools anaconda3/2.1.0 spades/3.13.1 perl/5.20.2 && metaspades.py -1 '+read1+' -2 '+read2+' -m '+memory+' -k '+k_spades+' --only-assembler -o '+out+'')
    subprocess.check_call(spadesCmd, shell=True)
