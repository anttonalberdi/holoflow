#13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-ia', help="index assembly file", dest="idx_a", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
idx_a=args.idx_a
ID=args.ID
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tAssembly Indexing step - '+ID+'\n')
    log.write('The assembly file needs to be indexed so the original read files can be mapped to it.\n\n')

# if the .fai indexed assembly file does not exist, continue
if not os.path.isfile(idx_a):

    # index assembly with samtools and bwa, both necessary for further steps
    idxsamCmd='module load tools samtools/1.11 && samtools faidx '+a+''
    idxbwaCmd='module load tools bowtie2/2.4.2 \
    && bowtie2-build --large-index --threads 40 '+a+' '+a+''

    subprocess.Popen(idxbwaCmd, shell=True).wait()
    subprocess.Popen(idxsamCmd, shell=True).wait()
