#13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-ia', help="index assembly file", dest="idx_a", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
idx_a=args.idx_a
sample=args.sample
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tAssembly Indexing step - Sample '+sample+'\n')
    log.write('The assembly file needs to be indexed so the original read files can be mapped to it.\n\n')


if not (os.path.exists(str(idx_a))):
    idxsamCmd='module load tools samtools/1.9 && samtools faidx '+a+''
    idxbwaCmd='module load tools bwa/0.7.15 && bwa index '+a+''

    subprocess.check_call(idxbwaCmd, shell=True)
    subprocess.check_call(idxsamCmd, shell=True)
