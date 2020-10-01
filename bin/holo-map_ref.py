#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-refg', help="reference genomes", dest="ref_gen", required=True)
parser.add_argument('-obam', help="all bam file", dest="all_bam", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-k', help="minimum seed length", dest="k", required=True)
parser.add_argument('-w', help="band width", dest="w", required=True)
parser.add_argument('-d', help="extension score threshold", dest="d", required=True)
parser.add_argument('-A', help="matching score", dest="A", required=True)
parser.add_argument('-B', help="mismatch penalty", dest="B", required=True)
parser.add_argument('-O', help="gap open penalty", dest="O", required=True)
parser.add_argument('-E', help="gap extension penalty", dest="E", required=True)
parser.add_argument('-L', help="clipping penalty", dest="L", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
#parser.add_argument('-R', help="Complete read group header line", dest="R", required=True)
args = parser.parse_args()

all_bam=args.all_bam
read1=args.read1
read2=args.read2
ref_gen=args.ref_gen
t=args.t
k=args.k
w=args.w
d=args.d
A=args.A
B=args.B
O=args.O
E=args.E
L=args.L
sample=args.sample
log=args.log
#R=args.R


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMapping To Reference Genomes step - Sample '+sample+'\n')
    log.write('All the reads are being mapped to the reference genome(s).\n')


if (k == "loose"):
    mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -k 19 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
    subprocess.check_call(mapCmd, shell=True)


if (k == "semistringent"):
    mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -k 30 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
    subprocess.check_call(mapCmd, shell=True)


if (k == "superstringent"):
    mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -k 50 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
    subprocess.check_call(mapCmd, shell=True)

if not ((k == "loose") or (k == "semistringent") or (k == "superstringent")):
    print(''+k+' is not a valid value, k = loose/semistringent/stringent - See config.yaml')
