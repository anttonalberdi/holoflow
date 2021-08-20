#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import gzip

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-refg', help="reference genomes", dest="ref_gen", required=True)
parser.add_argument('-ibam', help="all bam file", dest="all_bam", required=True)
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-obam', help="bam file", dest="bam", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-si', help="stats input file", dest="in_stats", required=True)
parser.add_argument('-so', help="stats output file", dest="out_stats", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
args = parser.parse_args()

all_bam=args.all_bam
ref_gen=args.ref_gen
bam=args.bam
read1=args.read1
read2=args.read2
log=args.log
in_stats=args.in_stats
out_stats=args.out_stats
ID=args.ID

# Run
# Write to log
with open(str(log),'a+') as logi:
    logi.write('A .bam file is generated containing the mapped reads, and two .fastq files containing the metagenomic ones.\n\n')

# sort bam for genomics
refbam1Cmd = 'module load tools samtools/1.11 && samtools view -T '+ref_gen+' -b -F12 '+all_bam+' > '+bam+'.notsorted && samtools sort -T '+bam+'.'+ID+' -o '+bam+' '+bam+'.notsorted && rm '+bam+'.notsorted'
subprocess.check_call(refbam1Cmd, shell=True)

# extract not-mapped to the reference genome reads + keep reference bam
refbam2Cmd = 'module load tools samtools/1.11 && samtools view -T '+ref_gen+' -b -f12 '+all_bam+' | samtools fastq -c 6 -1 '+read1+' -2 '+read2+' -'
subprocess.check_call(refbam2Cmd, shell=True)

# remove general bam
rmAllbamCmd = 'rm '+all_bam+'' # Change this if dark matter workflow
subprocess.check_call(rmAllbamCmd, shell=True)



    # Get stats after duplicate removal
mvstatsCmd= 'mv '+in_stats+' '+out_stats+''
subprocess.check_call(mvstatsCmd, shell=True)

reads = 0
bases = 0
with gzip.open(str(read1), 'rt') as read: # outputs are compressed files: .gz extension
    for id in read:
        seq = next(read)
        reads += 1
        bases += len(seq.strip())*2
        next(read)
        next(read)

#Print stats to statsfile
statsfile=open(str(out_stats),"a+")
statsfile.write("Reads after mapping to reference genome \t{0} ({1} bases)\r\n".format(reads,bases))
statsfile.close()


# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logo:
    logo.write('\t\t'+current_time+'\tPreprocessing with Holoflow has finished.\n')
