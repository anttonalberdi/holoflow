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
parser.add_argument('-M', help="picard-friendly bam", dest="picard", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
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
picard=args.picard
ID=args.ID
log=args.log
#R=args.R


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMapping To Reference Genomes step - '+ID+'\n')
    log.write('All the reads are being mapped to the reference genome(s).\n')

#de- compress inputs
if (os.path.exists(read1)):
    compressCmd1='gunzip '+read1+' '+read2+''
    subprocess.Popen(compressCmd1,shell=True).wait()
    read1 = read1.replace('.gz','')
    read2 = read2.replace('.gz','')

# not very optimal 
if (k == "loose"): # -k 19
    if not (picard == 'False'):
        mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -M -t '+t+' -k 19 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
        subprocess.check_call(mapCmd, shell=True)
    else:
        mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -t '+t+' -k 19 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
        subprocess.check_call(mapCmd, shell=True)


if (k == "semistringent"): # -k 30
    if not (picard == 'False'):
        mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -M -t '+t+' -k 30 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
        subprocess.check_call(mapCmd, shell=True)
    else:
        mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -t '+t+' -k 30 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
        subprocess.check_call(mapCmd, shell=True)


if (k == "superstringent"): # -k 50
    if not (picard == 'False'):
        mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -M -t '+t+' -k 50 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
        subprocess.check_call(mapCmd, shell=True)
    else:
        mapCmd = 'module load tools samtools/1.11 bwa/0.7.15 && bwa mem -t '+t+' -k 50 -w '+w+' -d '+d+' -A '+A+' -B '+B+' -O '+O+' -E '+E+' -L '+L+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+ref_gen+' '+read1+' '+read2+' | samtools view -T '+ref_gen+' -b - > '+all_bam+''
        subprocess.check_call(mapCmd, shell=True)

if not ((k == "loose") or (k == "semistringent") or (k == "superstringent")):
    print(''+k+' is not a valid value, k = loose/semistringent/stringent - See config.yaml')

# re -compress inputs
if (os.path.isfile(all_bam)):
    compressCmd2='gzip '+read1+' '+read2+''
    subprocess.Popen(compressCmd2,shell=True).wait()
