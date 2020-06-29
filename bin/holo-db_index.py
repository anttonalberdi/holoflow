#19.06.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-db', help="data base file", dest="db", required=True)
parser.add_argument('-idx_bwa', help="index data base file bwa", dest="idx_bwa", required=True)
parser.add_argument('-idx_smt', help="index data base file samtools", dest="idx_smt", required=True)
args = parser.parse_args()


db=args.db
idx_bwa=args.idx_bwa
idx_smt=args.idx_smt


# Run
if (os.path.exists(str(idx_bwa)) and os.path.exists(str(idx_smt))):
    pass

else:
    # first decompress db
    if str(db).endswith(".gz"):
        decompressCmd=('gunzip '+db+'')
        subprocess.check_call(decompressCmd, shell=True)
        decomp_db= db.replace('.gz','')

        # index
        idxsamCmd='module load tools samtools/1.9 && samtools faidx '+decomp_db+''
        idxbwaCmd='module load bwa/0.7.15 && bwa index '+decomp_db+''
        subprocess.check_call(idxbwaCmd, shell=True)
        subprocess.check_call(idxsamCmd, shell=True)

    else:
        # index
        idxsamCmd='module load tools samtools/1.9 && samtools faidx '+db+''
        idxbwaCmd='module load bwa/0.7.15 && bwa index '+db+''
        subprocess.check_call(idxbwaCmd, shell=True)
        subprocess.check_call(idxsamCmd, shell=True)
