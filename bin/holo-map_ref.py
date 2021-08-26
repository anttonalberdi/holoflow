#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-refg', help="reference genomes", dest="ref_gen", required=True)
parser.add_argument('-obam', help="all bam file", dest="all_bam", required=True)
parser.add_argument('-threads_bt2', help="threads to use", dest="threads", required=True)
# parser.add_argument('-sensitivity', help="Bowtie2 sensitivity setting, default sensitive", dest="sens", required=True)
# parser.add_argument('-alignmenttype', help="Bowtie2 alignment type, default end-to-end", dest="at", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

all_bam=args.all_bam
read1=args.read1
read2=args.read2
ref_gen=args.ref_gen
#picard=args.picard
threads=args.threads
# sens=args.sens
# at=args.at
ID=args.ID
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMapping To Reference Genomes step - '+ID+'\n')
    log.write('All the reads are being mapped to the reference genome(s).\n')


# if (k == "loose"): # -k 19
#     if not (picard == 'False'):
mapCmd = 'module load tools samtools/1.11 bowtie2/2.4.2 pigz/2.3.4 \
          && bowtie2 \
          --time \
          --threads '+threads+' \
          --rg-id "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:'+ID+'" \
          -x '+ref_gen+' \
          -1 '+read1+' \
          -2 '+read2+' | samtools view -@ '+threads+' -T '+ref_gen+' -b -o '+all_bam+''
subprocess.check_call(mapCmd, shell=True)
#     else:
#         mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -t '+t+' -k 19 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:'+ID+'" '+ref_gen+' <(gunzip -c '+read1+') <(gunzip -c '+read2+') | samtools view -T '+ref_gen+' -b - > '+all_bam+''
#         subprocess.check_call(mapCmd, shell=True,executable="/bin/bash")

# re -compress inputs
# if (os.path.isfile(all_bam)):
#     compressCmd2='pigz -p '+threads+' '+read1+' & pigz -p '+threads+' '+read2+''
#     subprocess.Popen(compressCmd2,shell=True,executable="/bin/bash").wait()
