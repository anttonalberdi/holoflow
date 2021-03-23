#28.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-o', help="output directory", dest="out", required=True)
parser.add_argument('-empty_o', help="empty touched file", dest="empty_o", required=True)
parser.add_argument('-coa', help='coassembly', dest="coassembly", required=False)
parser.add_argument('-m', help="memory", dest="memory", required=False)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-k_megahit', help="k-mer size list megahit", dest="k_megahit", required=True)
parser.add_argument('-k_spades', help="k-mer size list spades", dest="k_spades", required=False)
parser.add_argument('-a', help="assembler", dest="assembler", required=False)
parser.add_argument('-temp_a', help="temporal assembly file", dest="temp_a", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


read1=args.read1
read2=args.read2
out=args.out
k_megahit=args.k_megahit
threads=args.threads
empty_o=args.empty_o
temp_a=args.temp_a
ID=args.ID
log=args.log


if (args.coassembly):
    args.assembler='megahit'
    assembler=args.assembler

# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\tHOLOFLOW\tMETAGENOMICS\n\t\t'+current_time+'\tMetagenomic Data Assembly step - '+ID+'\n')
    log.write('The .fastq files coming from Holoflow Preprocessing, are those which could not be mapped to a \nreference genome. These contain the metagenomic reads; as no reference genome exists to them,\n they have to be assembled de novo. This is done by '+args.assembler+' here, which sorts the reads together into\ncontigs or scaffolds giving out one only assembly fasta file.\n\n')


if os.path.exists(temp_a):
    pass

if not os.path.exists(temp_a):

    if (args.assembler == "megahit"):

        if (args.coassembly):

            with open(read1,'r') as f1, open(read2,'r') as f2:
                read1_paths = f1.readline()
                read2_paths = f2.readline()

            megahitCmd = 'module load tools megahit/1.2.9 && megahit -1 '+read1_paths+' -2 '+read2_paths+' -t '+threads+' --k-list '+k_megahit+' -o '+out+''
            subprocess.check_call(megahitCmd, shell=True)

            mv_megahitCmd = 'mv '+out+'/final.contigs.fa '+out+'/temp_assembly.fa'
            subprocess.check_call(mv_megahitCmd, shell=True)

        else:

            megahitCmd = 'module load tools megahit/1.2.9 && megahit -1 '+read1+' -2 '+read2+' -t '+threads+' --k-list '+k_megahit+' -o '+out+''
            subprocess.check_call(megahitCmd, shell=True)

            mv_megahitCmd = 'mv '+out+'/final.contigs.fa '+out+'/temp_assembly.fa'
            subprocess.check_call(mv_megahitCmd, shell=True)


    if args.assembler == "spades":

        if (args.coassembly):

            with open(read1,'r') as f1, open(read2,'r') as f2:
                read1_paths = f1.readline()
                read2_paths = f2.readline()

                # Merge all read1, read2's content into 1 file each
                read1_coa = read1.replace('_1.fastq','merged_1.fastq')
                read2_coa = read1.replace('_2.fastq','merged_2.fastq')

                mergeCmd = 'cat '+read1_paths+' > '+read1_coa+' && cat '+read2_paths+' > '+read2_coa+''
                subprocess.check_call(mergeCmd, shell=True)

                # Run spades on merged files
                spadesCmd = 'module unload anaconda3/4.4.0 && mkdir '+out+' && module load tools anaconda3/2.1.0 spades/3.13.1 perl/5.20.2 && metaspades.py -1 '+read1_coa+' -2 '+read2_coa+' -m '+args.memory+' -k '+args.k_spades+' --only-assembler -o '+out+''
                subprocess.check_call(spadesCmd, shell=True)

                mv_spadesCmd = 'mv '+out+'/scaffolds.fasta '+out+'/temp_assembly.fa'
                subprocess.check_call(mv_spadesCmd, shell=True)


        else:

            spadesCmd = 'module unload anaconda3/4.4.0 && mkdir '+out+' && module load tools anaconda3/2.1.0 spades/3.13.1 perl/5.20.2 && metaspades.py -1 '+read1+' -2 '+read2+' -m '+args.memory+' -k '+args.k_spades+' --only-assembler -o '+out+''
            subprocess.check_call(spadesCmd, shell=True)

            mv_spadesCmd = 'mv '+out+'/scaffolds.fasta '+out+'/temp_assembly.fa'
            subprocess.check_call(mv_spadesCmd, shell=True)


    emptytouchCmd='touch '+empty_o+''
    subprocess.check_call(emptytouchCmd, shell=True)
