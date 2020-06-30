#13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-i', help="input assembly file", dest="i", required=True)
parser.add_argument('-o', help="output genetic coordinates", dest="o", required=True)
parser.add_argument('-a', help="protein translations", dest="a", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

i=args.i
o=args.o
a=args.a
sample=args.sample
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tProdigal Protein Prediction step - Sample '+sample+'\n')
    log.write('Prodigal is a gene-finding program for microbial sequences, which will be used in following taxonomic\nassignation procedures.\n\n')


prodigalCmd='module unload gcc && module load tools prodigal/2.6.3 && prodigal -i '+i+' -o '+o+' -a '+a+' -p meta'
subprocess.check_call(prodigalCmd, shell=True)
