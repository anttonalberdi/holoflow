#29.04.2021 - Holoflow

import subprocess
import argparse
import os
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-out_dir', help="output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
out_dir=args.out_dir
ID=args.ID
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tAssembly Annotation step - '+ID+'\n')
    log.write('The assembly file is being functionally annotated by DRAM v1.2.0 (Distilled and Refined Annotation of Metabolism).\nFirst an annotation step to assign database identifiers to gene and then a distill step to curate these annotations into useful functional categories.\n\n')

# Run annotation
if os.path.isfile(a):
    dram1Cmd='module load dram/1.2.0 && DRAM.py annotate -i '+a+' -o '+out_dir+''
    subprocess.Popen(dram1Cmd,shell=True).wait()

# In the output annotation folder there will be various files. genes.faa and genes.fna are fasta files with all genes called by prodigal
# with additional header information gained from the annotation as nucleotide and amino acid records respectively. genes.gff is a GFF3
# with the same annotation information as well as gene locations. scaffolds.fna is a collection of all scaffolds/contigs given as input
# to DRAM.py annotate with added bin information in the headers. annotations.tsv is the most important output of the annotation. This
# includes all annotation information about every gene from all MAGs. Each line is a different gene and each column contains annotation
# information. trnas.tsv contains a summary of the tRNAs found in each MAG.

    # Summarise annotation
    dram2Cmd='DRAM.py distill -i '+out_dir+'/annotations.tsv -o '+out_dir+'/summary --trna_path '+out_dir+'/trnas.tsv --rrna_path '+out_dir+'/rrnas.tsv'
    #subprocess.Popen(dram1Cmd,shell=True).wait()
