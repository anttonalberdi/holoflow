#19.06.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-db', help="data base file", dest="db", required=True)
parser.add_argument('-check', help="file OK", dest="check", required=True)
args = parser.parse_args()


db=args.db
check=args.check


#   d. If all files are FINE, create tiny .txt which says it worked, just x checking, if not:BREAK
#   e. Compress ALL OUTPUT FILES outputdir (-d)/NameOutputDB.fna.tar.gz

# Run
if (os.path.exists(str(idx_db))): # if fasta has been correctly assembled
    file = os.path.dirname(sys.argv[0])
    curr_dir = os.path.abspath(file)

    compressCmd=('tar -zcvf '+db+'.tar.gz '+curr_dir+'')
    subprocess.check_call(compressCmd, shell=True)

    with open(str(check),'w') as check_file:
        check_file.write('All reference genomes have been merged and indexed successfully.')
