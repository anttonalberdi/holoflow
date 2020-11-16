#16.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import os
import sys


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bt_mtb', help="metabat bin table", dest="bt_mtb", required=True)
parser.add_argument('-bt_mxb', help="maxbin bin table", dest="bt_mxb", required=True)
parser.add_argument('-mtb', help="metabat bin dir", dest="mtb", required=True)
parser.add_argument('-mxb', help="maxbin bin dir", dest="mxb", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


bt_mtb=args.bt_mtb
bt_mxb=args.bt_mxb
mtb=args.mtb
mxb=args.mxb
ID=args.ID
log=args.log

##############################################
#################### WRITE TO LOG ##########################
##############################################

    # If only one of the binners produced bins:

bt_todupl=''
bp_todupl=''
bt_e=''
bp_e=''
dupl_binner=''
empty_binner=''

if not (os.path.exists(bt_mtb)):
    bt_todupl=bt_mxb
    bp_todupl=mxb
    dupl_binner='mxb'

    bt_e=bt_mtb
    bp_e=mtb
    empty_binner='mtb'

if not (os.path.exists(bt_mxb)):
    bt_todupl=bt_mtb
    bp_todupl=mtb
    dupl_binner='mtb'

    bt_e=bt_mxb
    bp_e=mxb
    empty_binner='mxb'

if (os.path.exists(bt_mxb) and os.path.exists(bt_mtb)):
    sys.exit()


# Duplicate the existing bins and bin table and rename duplicates
if os.path.exists(bp_e):
    os.rmdir(bp_e)

    mvCmd='cp -r '+bp_todupl+' '+bp_e+' && cp '+bt_todupl+' '+bt_e+' && grep '+str(dupl_binner)+' '+bp_e+'/* | for f in ; do mv "$f" "$(echo "$f" | sed s/'+str(dupl_binner)+'/'+str(empty_binner)+'/)"; done && grep '+str(dupl_binner)+' '+bt_e+' | sed s/'+str(dupl_binner)+'/'+str(empty_binner)+'/'
    subprocess.check_call(mvCmd,shell=True)

else:
    pass



    
