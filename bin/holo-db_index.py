#19.06.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-db', help="data base file", dest="db", required=True)
parser.add_argument('-idx_bwa', help="index data base file bwa", dest="idx_bwa", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-idx_smt', help="index data base file samtools", dest="idx_smt", required=True)
args = parser.parse_args()


db=args.db
log=args.log
idx_bwa=args.idx_bwa
idx_smt=args.idx_smt


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),"w+") as log:
    log.write('\tHOLOFLOW\tPREPARE GENOMES\n\t\t'+current_time+'\tData Base indexing step\n')
    log.write('The data base needs to be indexed with BWA and SAMTOOLS so the reads can be mapped to it\nduring preprocessing.\n\n')



# first decompress db if necessary
if str(db).endswith(".gz"):
    decompressCmd=('gunzip '+db+'')
    subprocess.check_call(decompressCmd, shell=True)
    decomp_db= db.replace('.gz','')

else:
    decomp_db = db

# Index

if os.path.exists(str(idx_bwa)):
    pass

else:
    idxbwaCmd='module load tools bwa/0.7.15 && bwa index '+decomp_db+'' ###### bwa cores 1 
    subprocess.check_call(idxbwaCmd, shell=True)



if os.path.exists(str(idx_smt)):
    pass

else:
    # index
    idxsamCmd='module load tools samtools/1.11 && samtools faidx '+decomp_db+''
    subprocess.check_call(idxsamCmd, shell=True)
