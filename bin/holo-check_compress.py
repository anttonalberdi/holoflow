#19.06.2020 - Holoflow 0.1.

import subprocess
import argparse
import sys
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-db', help="data base file", dest="db", required=True)
parser.add_argument('-idx_db', help="indexed data base file", dest="idx_db", required=True)
parser.add_argument('-db_dir', help="data base directory", dest="db_dir", required=True)
parser.add_argument('-check', help="file OK", dest="check", required=True)
args = parser.parse_args()


db=args.db
idx_db=args.idx_db
db_dir=args.db_dir
check=args.check



# Run
if (os.path.exists(str(idx_db)) and os.path.exists(str(db))):

    compressCmd=('tar -zcvf '+db+'.tar.gz '+db_dir+'')
    subprocess.check_call(compressCmd, shell=True)

    with open(str(check),'w') as check_file:
        check_file.write('All reference genomes have been merged and indexed successfully.')
