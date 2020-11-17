#16.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import os
import sys


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-binning_dir', help="binning directory", dest="binning_dir", required=True)
parser.add_argument('-check_file', help="empty check file", dest="check_file", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


binning_dir=args.binning_dir
check_file=args.check_file
ID=args.ID
log=args.log

##############################################
#################### WRITE TO LOG ##########################
##############################################

mtb=str(os.path.join(binning_dir,ID+'_metabat'))
bt_mtb=str(binning_dir+'/'+ID+'.bins_metabat.txt')
mxb=str(os.path.join(binning_dir,ID+'_maxbin'))
bt_mxb=str(binning_dir+'/'+ID+'.bins_maxbin.txt')



    # If only one of the binners produced bins:
bt_todupl=''
bp_todupl=''
bt_e=''
bp_e=''
dupl_binner=''
empty_binner=''

if not (os.path.isfile(bt_mtb)):
    bt_todupl=bt_mxb
    bp_todupl=mxb
    dupl_binner='mxb'

    bt_e=bt_mtb
    bp_e=mtb
    empty_binner='mtb'

    if os.path.exists(bp_e):
        os.rmdir(bp_e)

if not (os.path.isfile(bt_mxb)):
    bt_todupl=bt_mtb
    bp_todupl=mtb
    dupl_binner='mtb'

    bt_e=bt_mxb
    bp_e=mxb
    empty_binner='mxb'

    if os.path.exists(bp_e):
        os.rmdir(bp_e)

else:
    os.mknod(str(check_file))
    sys.exit()


# Duplicate the existing bins and bin table and rename duplicates
mvCmd='cp -r '+bp_todupl+' '+bp_e+' && for f in '+bp_e+'/*'+str(dupl_binner)+'* ; do mv "$f" "$(echo "$f" | sed s/'+str(dupl_binner)+'/dup_'+str(empty_binner)+'/)"; done'
subprocess.Popen(mvCmd,shell=True).wait()

cpCmd='cp '+bt_todupl+' '+bt_e+'.tmp && grep '+str(dupl_binner)+' '+bt_e+'.tmp | sed s/'+str(dupl_binner)+'/dup_'+str(empty_binner)+'/ > '+bt_e+' && rm '+bt_e+'.tmp'
subprocess.Popen(cpCmd,shell=True).wait()

emptyCmd='touch '+check_file+''
subprocess.Popen(emptyCmd,shell=True).wait()
