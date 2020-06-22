#19.06.2020 - Holoflow 0.1.

import argparse
import subprocess
import glob
import os
import sys
import ruamel.yaml


###########################
#Argument parsing
###########################
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-refp', help="path to reference genomes", dest="refp", required=True)
parser.add_argument('-suff', help="reference genomes common termination", dest="suff", required=True)
parser.add_argument('-DB', help="data base file name", dest="DB", required=True)
args = parser.parse_args()

refp=args.refp
suff=args.suff
DB=args.DB


# obtain full paths of files
ref_genomes = os.path.abspath(x) for x in glob.glob(''+refp+'/*'+suff+'')

# reformat genomes

# merge genomes
mergeCmd=(''+ref_genomes+' > '+DB+'')
subprocess.check_call(mergeCmd, shell=True)
