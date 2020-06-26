#19.06.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-db', help="data base file", dest="db", required=True)
parser.add_argument('-idb', help="index data base file", dest="idx_db", required=True)
args = parser.parse_args()


db=args.db
idx_db=args.idx_db


# Run
if not (os.path.exists(str(idx_db))):
    # first decompress db
    if str(db).endswith(".gz"):
        decompressCmd=('gunzip '+db+'')
        subprocess.check_call(decompressCmd, shell=True)
        decomp_db= db.replace('.gz','')

    else:
        decomp_db = db

    # index
    idxsamCmd='module load tools samtools/1.9 && samtools faidx '+decomp_db+''
    idxbwaCmd='module load bwa/0.7.15 && bwa index '+decomp_db+''

    subprocess.check_call(idxbwaCmd, shell=True)
    subprocess.check_call(idxsamCmd, shell=True)
    
else:
    pass
