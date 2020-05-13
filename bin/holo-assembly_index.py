#13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-ia', help="index assembly file", dest="idx_a", required=True)
args = parser.parse_args()


a=args.a
idx_a=args.idx_a


# Run
if not os.path.exists(str(idx_a)):
    idxCmd='module load tools samtools/1.9 && samtools faidx '+a+' && module load tools bwa/0.7.15 && bwa index '+a+''
    subprocess.check_call(idxCmd, shell=True)
else:
    pass
