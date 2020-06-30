#14.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-mtb', help="metabat depth file", dest="mtb", required=True)
parser.add_argument('-mxb', help="maxbin depth file", dest="mxb", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
mtb=args.mtb
mxb=args.mxb
sample=args.sample
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tDepth File Generation step - Sample '+sample+'\n')
    log.write('Depth file containing coverage info about the reads is being generated to be used during binning.\n\n')


# Metabat
metabatCmd='module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+mtb+' '+a+''
subprocess.check_call(metabatCmd, shell=True)


# Maxbin
maxbinCmd='cp '+mtb+' '+mxb+''
subprocess.check_call(maxbinCmd, shell=True)
