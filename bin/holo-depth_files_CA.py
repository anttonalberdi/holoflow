#14.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-mtb', help="metabat depth file", dest="mtb", required=True)
parser.add_argument('-mxb', help="maxbin depth file", dest="mxb", required=True)
parser.add_argument('-cct', help="concoct depth file", dest="cct", required=True)
args = parser.parse_args()


a=args.a
mtb=args.mtb
mxb=args.mxb
cct=args.cct


# Run

# Metabat
metabatCmd='module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+mtb+' '+a+''
subprocess.check_call(metabatCmd, shell=True)


# Maxbin
maxbinCmd='cp '+mtb+' '+mxb+''
subprocess.check_call(maxbinCmd, shell=True)

#Concoct
concoctCmd='cat '+mtb+' | cut -f-1,4,6,8- > '+cct+''
subprocess.check_call(concoctCmd, shell=True)
