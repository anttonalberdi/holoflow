#19.06.2020 - Holoflow 0.1.

import subprocess
import argparse
import sys
import time
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-db', help="data base path", dest="db", required=True)
parser.add_argument('-idx_db', help="indexed data base file", dest="idx_db", required=True)
parser.add_argument('-db_dir', help="data base directory", dest="db_dir", required=True)
parser.add_argument('-db_ID', help="data base ID", dest="db_ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-check', help="file OK", dest="check", required=True)
args = parser.parse_args()


db=args.db
idx_db=args.idx_db
db_dir=args.db_dir
db_ID=args.db_ID
log=args.log
check=args.check



# Run
# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\t\t'+current_time+'\tCompressing data base and index files step\n\n')
    logi.close()


if (os.path.exists(str(idx_db)) and os.path.exists(str(db))) and (not os.path.exists(str(check))):

    with open(str(check),'w') as check_file:
        check_file.write('All reference genomes have been merged and indexed successfully.')

    compressCmd=('cd '+db_dir+' && tar -zcvf ../'+db_ID+'.tar.gz * && rm -rf '+db_dir+'')
    subprocess.check_call(compressCmd, shell=True)


# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logf:
    logf.write('\t\t'+current_time+'\tHoloflow has completed the preparation of the reference genomes.\n\n')
