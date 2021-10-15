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
parser.add_argument('-db_names', help="names of the db/dbs to be used", dest="db_names", required=True)
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
    logi.write('Genes which map to the designed database(s) {Plants, Invertebrates...} will be annotated by Diamond 2.0.6.\n\n')

if not os.path.exists(out_dir):
    mkdirCmd='mkdir -p '+out_dir+''
    subprocess.Popen(mkdirCmd,shell=True).wait()

            ####################
            #### MERGED dbs option ---> COMPETITIVE
            ####################
# merge all db that the user wants to map the predicted ORFs to
tmp_dbs = out_dir+'/'+db_names+'-TMP_merge.dat.gz'

if not os.path.isfile(tmp_dbs):
        # find dbs in db dir
    db_files = glob.glob(db_dir+'/*.fasta.gz')
    db_tomerge = ''
     # generate a string with those dbs to merge
    for db_path in db_files:
        # find all databases in db dir
        for db_name in db_names.strip().split('_'):
            # get names of the tax. groups the user wants to annotate from, _ delim
            if db_name in db_path:
                db_tomerge += db_path+' '       # create string with paths to selected dbs

            else:
                pass
    mergeCmd='cat '+db_tomerge+' > '+tmp_dbs+''    # merge the selected dbs into one file
    subprocess.Popen(mergeCmd,shell=True).wait()


# annot
if os.path.isfile(tmp_dbs):
    out_annot = out_dir+'/'+db_names+'-annotation.dmnd'

    diamondCmd='module load diamond/2.0.6 && \
    diamond blastp -d '+tmp_dbs+' -q '+faa+' -o '+out_annot+' -p '+t+' -k 1'
    subprocess.Popen(diamondCmd, shell=True).wait()
    # given the database and the predicted proteins, retrieve the first best map match


#             ####################
#             #### ONE DB BY ONE
#             ####################
#
# # find dbs in db dir
# db_files = glob.glob(db_dir+'/*.fasta.gz')
# db_toannot = list()
#  # generate a string with those dbs to merge
# for db_path in db_files:                    # find all databases in db dir
#     for db_name in db_names.split('_'):     # get names of the tax. groups the user wants to annotate from, _ delim
#         if db_name in db_path:
#             db_toannot.append(db_path.strip())      # create list with dbs paths
#         else:
#             pass
#
# # annot
# for db_annot in db_toannot:
#     db_name = db_annot.replace(db_dir,'').replace('.fasta.gz','')
#
#     if os.path.isfile(db_annot):
#         out_annot = out_dir+'/'+ID+'-'+db_name+'_annot.dmnd'
#
#         diamondCmd='module load diamond/2.0.6 && diamond blastp -d '+db_annot+' -q '+faa+' -o '+out_annot+' -p '+t+' -k 1'
#         subprocess.Popen(diamondCmd, shell=True).wait()
