#09.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-megahit', help="megahit input", dest="megahit", required=True)
parser.add_argument('-spades', help="spades input", dest="spades", required=True)
parser.add_argument('-o', help="output directory", dest="out", required=True)
parser.add_argument('-a', help="assembler", dest="assembler", required=True)
args = parser.parse_args()


out=args.out
megahit=args.megahit
spades=args.spades
assembler=args.assembler

# Run
if assembler == "megahit":
    megahitCmd = shell('mv '+megahit+' '+out+'')
    subprocess.check_call(megahitCmd, shell=True)

if assembler == "spades":
    spadesCmd = shell('mv '+spades+' '+out+'')
    subprocess.check_call(spadesCmd, shell=True)
