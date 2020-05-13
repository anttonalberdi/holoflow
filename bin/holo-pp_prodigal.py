#13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-i', help="input assembly file", dest="a", required=True)
parser.add_argument('-o', help="output genetic coordinates", dest="o", required=True)
parser.add_argument('-a', help="protein translations", dest="a", required=True)
args = parser.parse_args()

i=args.i
o=args.o
a=args.a


# Run
prodigalCmd='module unload gcc && module load tools prodigal/2.6.3 && prodigal -i '+i+' -o '+o+' -a '+a+' -p meta'
subprocess.check_call(prodigalCmd, shell=True)
