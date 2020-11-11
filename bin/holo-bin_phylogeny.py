#11.11.2020

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-gen_dir', help="prokka genes directory", dest="gene_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


gene_dir=args.gene_dir
out_dir=args.out_dir
ID=args.ID
log=args.log
threads=args.threads


# Run
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tTaxonomic Classification step - '+ID+'\n')
        logi.write('GTDB-Tk is assigning objective taxonomic classifications to baterial genomes based on the Genome Database Taxonomy GTDB.\n\n')

        
