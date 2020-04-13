#08.04.2020 - Holoflow 0.1

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-i', help="input_all", dest="input", required=True)
parser.add_argument('-sep', help="sep", dest="separator", required=True)
args = parser.parse_args()

input_file=args.input
read1=args.read1
read2=args.read2
separator=args.separator

# Run
cut1Cmd = 'cut --delimiter='+str(separator)+' -f1 '+input_file+' > '+read1+''
subprocess.check_call(cut1Cmd, shell=True)
cut2Cmd = 'cut --delimiter='+str(separator)+' -f2 '+input_file+' > '+read2+''
subprocess.check_call(cut2Cmd, shell=True)
rmCmd = 'rm '+input_file+''
subprocess.check_call(rmCmd, shell=True)
