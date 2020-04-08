#08.04.2020 - Holoflow 0.1

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-o1', help="path1", dest="read1", required=True)
parser.add_argument('-o2', help="path2", dest="read2", required=True)
parser.add_argument('-i', help="input_all", dest="input", required=True)
parser.add_argument('-sep', help="sep", dest="separator", required=True)
args = parser.parse_args()

input_file=args.input
read1=args.read1
read2=args.read2
separator=args.separator

# Run
cutCmd = 'cut --delimiter=separator -f1'+input_file+' > '+read1+' && cut --delimiter='+separator+' -f2 '+input+' > '+read2+' && rm '+input'
subprocess.check_call(cutCmd, shell=True)
