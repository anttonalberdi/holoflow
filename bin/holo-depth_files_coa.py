#03.12.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bam_p', help="path to bam files", dest="bam_p", required=True)
parser.add_argument('-mtb', help="metabat depth file", dest="mtb", required=True)
parser.add_argument('-mxb', help="maxbin depth file", dest="mxb", required=True)
parser.add_argument('-cct', help="concoct depth file ", dest="cct", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


bam_p=args.bam_p
mtb=args.mtb
mxb=args.mxb
cct=args.cct
ID=args.ID
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tDepth File Generation step - '+ID+'\n')
    log.write('Depth file containing coverage info about the reads is being generated to be used during binning.\n\n')



# Metabat
if not (os.path.isfile(mtb)):
    metabatCmd='module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+mtb+' '+bam_p+'/*.bam'
    subprocess.check_call(metabatCmd, shell=True)

# Concoct
if not (os.path.isfile(cct)):
    concoctCmd='cat '+mtb+' | cut -f1,4,6 > '+cct+''
    subprocess.Popen(concoctCmd, shell=True).wait()


# Maxbin
maxbinCmd='cp '+mtb+' '+mxb+''
subprocess.check_call(maxbinCmd, shell=True)
