#06.05.2021 - Holoflow 0.1.
import subprocess
import argparse
import os
import time
import glob


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-faa', help="protein sequences predicted ORFs", dest="faa", required=True)
parser.add_argument('-db_dir', help="db directory", dest="db_dir", required=True)
parser.add_argument('-db_names', help="names of the db/dbs to be used", dest="out_dir", required=True)
parser.add_argument('-out_dir', help="out_dir", dest="out_dir", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


faa=args.faa
db_names=args.db_names
db_dir=args.db_dir
out_dir=args.out_dir
t=args.threads
ID=args.ID
log=args.log



# Run
# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\tHOLOFLOW\tMETAGENOMICS\n\t\t'+current_time+'\t - '+ID+'\n')
    logi.write('   \n\n')

# merge all db that the user wants to map the predicted ORFs to
tmp_dbs = out_dir+'/'+db_names+'-TMP_merge.dat.gz' # don't know if it is better to merge them or would be better to 1 by 1

if not os.path.isfile(tmp_dbs):
        # find dbs in db dir
    db_files = glob.glob(db_dir+'/*.dat.gz')
    db_tomerge = ''
     # generate a string with those dbs to merge
    for db_path in db_files:
        for db_name in db_names.split('_'):
            if db_name in db_path:
                db_tomerge += db_path+' '
            else:
                pass

    mergeCmd='zcat '+db_tomerge+' > tmp_dbs'
    subprocess.Popen(mergeCmd,shell=True).wait()

# annot
if os.path.isfile(tmp_dbs):
    out_annot = out_dir+'/'+db_names+'-annotation.dmnd'

    diamondCmd='module load diamond/2.0.6 && diamond blastp -d '+tmp_dbs+' -q '+faa+' -o '+out_annot+' -p '+t+' -k 1'
    subprocess.Popen(diamondCmd, shell=True).wait()
