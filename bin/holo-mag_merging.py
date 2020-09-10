#10.09.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-dt_bd', help="dastool bin directory", dest="dt_bd", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()



dt_bd=args.dt_bd
out_dir=args.out_dir
sample=args.sample
log=args.log
threads=args.threads


# Run

sspaceDepCmd='module load tools perl/5.24.0 sspace-standard/3.0 parallel/20190522'
subprocess.check_call(sspaceDepCmd, shell=True)
