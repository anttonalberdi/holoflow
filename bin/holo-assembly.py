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


# if (args.coassembly):
#     args.assembler='megahit'
#     assembler=args.assembler

# Run
# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\tHOLOFLOW\tMETAGENOMICS\n\t\t'+current_time+'\tMetagenomic Data Assembly step - '+ID+'\n')
    logi.write('The .fastq files coming from Holoflow Preprocessing, are those which could not be mapped to a \nreference genome. These contain the metagenomic reads; as no reference genome exists to them,\n they have to be assembled de novo. This is done by '+args.assembler+' here, which sorts the reads together into\ncontigs or scaffolds giving out one only assembly fasta file.\n\n')


if os.path.exists(temp_a):
    pass

if not os.path.exists(temp_a):

    if (args.assembler == "megahit"): # MEGAHIT is OK with compressed input

        if (args.coassembly):

            with open(read1,'r') as f1, open(read2,'r') as f2:
                read1_paths = f1.readline()
                read2_paths = f2.readline()

            megahitCmd = 'module load tools megahit/1.2.9 && megahit -1 '+read1_paths+' -2 '+read2_paths+' -t '+threads+' --k-list '+k_megahit+' -o '+out+''
            subprocess.Popen(megahitCmd, shell=True).wait()

            mv_megahitCmd = 'mv '+out+'/final.contigs.fa '+out+'/temp_assembly.fa'
            subprocess.Popen(mv_megahitCmd, shell=True).wait()

        else:

            megahitCmd = 'module load tools megahit/1.2.9 && megahit -1 '+read1+' -2 '+read2+' -t '+threads+' --k-list '+k_megahit+' -o '+out+''
            subprocess.Popen(megahitCmd, shell=True).wait()

            mv_megahitCmd = 'mv '+out+'/final.contigs.fa '+out+'/temp_assembly.fa'
            subprocess.Popen(mv_megahitCmd, shell=True).wait()


    if args.assembler == "spades":

        if not os.path.exists(out):
            os.makedirs(out)

        if (args.coassembly):

            with open(read1,'r') as f1, open(read2,'r') as f2:
                read1_paths = f1.readline().strip().split(',')
                read1_paths = (' ').join(read1_paths)
                read2_paths = f2.readline().strip().split(',')
                read2_paths = (' ').join(read2_paths)

            # Merge all read1, read2's content into 1 file each
            if '.gz' in read1_paths:
                read1_coa = out+'/'+ID+'.merged_1.fastq.gz'
                read2_coa = out+'/'+ID+'.merged_2.fastq.gz'

                if not os.path.isfile(read1_coa):
                    mergeCmd = 'zcat '+read1_paths+' > '+read1_coa+' && zcat '+read2_paths+' > '+read2_coa+''
                    subprocess.Popen(mergeCmd, shell=True).wait()

            else:
                read1_coa = out+'/'+ID+'.merged_1.fastq'
                read2_coa = out+'/'+ID+'.merged_2.fastq'

                if not os.path.isfile(read1_coa):
                    mergeCmd = 'cat '+read1_paths+' > '+read1_coa+' && cat '+read2_paths+' > '+read2_coa+''
                    subprocess.Popen(mergeCmd, shell=True).wait()

            # Run spades on merged files
            spadesCmd = 'module unload anaconda3/4.4.0 && module load tools anaconda3/2.1.0 spades/3.13.1 perl/5.20.2 && metaspades.py -1 '+read1_coa+' -2 '+read2_coa+' -m '+args.memory+' -k '+args.k_spades+' --only-assembler -o '+out+''
            subprocess.Popen(spadesCmd, shell=True).wait()

            mv_spadesCmd = 'mv '+out+'/scaffolds.fasta '+out+'/temp_assembly.fa'
            subprocess.Popen(mv_spadesCmd, shell=True).wait()


        else:

            spadesCmd = 'module unload anaconda3/4.4.0 && module load tools anaconda3/2.1.0 spades/3.13.1 perl/5.20.2 && metaspades.py -1 '+read1+' -2 '+read2+' -m '+args.memory+' -k '+args.k_spades+' --only-assembler -o '+out+''
            subprocess.Popen(spadesCmd, shell=True).wait()

            mv_spadesCmd = 'mv '+out+'/scaffolds.fasta '+out+'/temp_assembly.fa'
            subprocess.Popen(mv_spadesCmd, shell=True).wait()


    emptytouchCmd='touch '+empty_o+''
    subprocess.Popen(emptytouchCmd, shell=True).wait()
