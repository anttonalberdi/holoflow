#29.04.2021 - Holoflow

import subprocess
import argparse
import os
import time
import sys


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-conda_env_file', help="conda_env_file", dest="conda_env_file", required=True)
parser.add_argument('-config', help="config to load dbs", dest="config_dbs", required=True)
parser.add_argument('-out_dir', help="output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-min_c_size', help="minimum contig size", dest="min_c_size", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
conda_env_file=args.conda_env_file
config_dbs=args.config_dbs
out_dir=args.out_dir
ID=args.ID
log=args.log
t=args.t
min_c_size=args.min_c_size


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tAssembly Annotation step - '+ID+'\n')
    log.write('The assembly file is being functionally annotated by DRAM v1.2.0 (Distilled and Refined Annotation of Metabolism).\nFirst an annotation step to assign database identifiers to gene and then a distill step to curate these annotations into useful functional categories.\n\n')

# Run annotation
if os.path.isfile(a):

# In the output annotation folder there will be various files. genes.faa and genes.fna are fasta files with all genes called by prodigal
# with additional header information gained from the annotation as nucleotide and amino acid records respectively. genes.gff is a GFF3
# with the same annotation information as well as gene locations. scaffolds.fna is a collection of all scaffolds/contigs given as input
# to DRAM.py annotate with added bin information in the headers. annotations.tsv is the most important output of the annotation. This
# includes all annotation information about every gene from all MAGs. Each line is a different gene and each column contains annotation
# information. trnas.tsv contains a summary of the tRNAs found in each MAG.

    # Call Rscript to generate sub-trees
    file = os.path.dirname(sys.argv[0])
    curr_dir = os.path.abspath(file)

    dram_conda_runCmd= 'bash '+curr_dir+'/holo-assembly_annotation.sh '+conda_env_file+' '+config_dbs+' '+a+' '+out_dir+' '+t+' '+min_c_size+''
    subprocess.Popen(dram_conda_runCmd,shell=True).wait()
