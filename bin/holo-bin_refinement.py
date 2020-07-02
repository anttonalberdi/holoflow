#01.07.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bin_dir', help="assembly file", dest="a", required=True)
parser.add_argument('-out_dir', help="metabat bin table", dest="bt_mtb", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

bin_dir=args.bin_dir
out_dir=args.out_dir
sample=args.sample
log=args.log



# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tRefineM Bin Refinement step - Sample '+sample+'\n')
    log.write('\n\n')
