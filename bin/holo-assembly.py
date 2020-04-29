#28.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-o', help="output directory", dest="out", required=True)
parser.add_argument('-empty_o', help="empty touched file", dest="empty_o", required=True)
parser.add_argument('-m', help="memory", dest="memory", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-k_megahit', help="k-mer size list megahit", dest="k_megahit", required=True)
parser.add_argument('-k_spades', help="k-mer size list spades", dest="k_spades", required=True)
parser.add_argument('-a', help="assembler", dest="assembler", required=True)
parser.add_argument('-temp_a', help="temporal assembly file", dest="temp_a", required=True)
args = parser.parse_args()


read1=args.read1
read2=args.read2
out=args.out
memory=args.memory
k_megahit=args.k_megahit
k_spades=args.k_spades
threads=args.threads
assembler=args.assembler
empty_o=args.empty_o
temp_a=args.temp_a


# Run
if not os.path.exists(str(out)):

    emptytouchCmd='touch '+empty_o+''
    subprocess.check_call(emptytouchCmd, shell=True)

    if assembler == "megahit":
        megahitCmd = 'module load tools megahit/1.1.1 && mkdir '+out+' && megahit -1 '+read1+' -2 '+read2+' -t '+threads+' --k-list '+k_megahit+' -o '+out+''
        subprocess.check_call(megahitCmd, shell=True)


        mv_megahitCmd = 'mv '+out+'/final.contigs.fa '+out+'/temp_assembly.fa'
        subprocess.check_call(mv_megahitCmd, shell=True)

    if assembler == "spades":
        spadesCmd = 'module unload anaconda3/4.4.0 && mkdir '+out+' && module load tools anaconda3/2.1.0 spades/3.13.1 perl/5.20.2 && metaspades.py -1 '+read1+' -2 '+read2+' -m '+memory+' -k '+k_spades+' --only-assembler -o '+out+''
        subprocess.check_call(spadesCmd, shell=True)

        mv_spadesCmd = 'mv '+out+'/scaffolds.fasta '+out+'/temp_assembly.fa'
        subprocess.check_call(mv_spadesCmd, shell=True)
