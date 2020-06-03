#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import gzip
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-i1', help="path1 input", dest="read1i", required=True)
parser.add_argument('-i2', help="path2 input", dest="read2i", required=True)
parser.add_argument('-o1', help="path1 output", dest="read1o", required=True)
parser.add_argument('-o2', help="path2 output", dest="read2o", required=True)
parser.add_argument('-a1', help="adapter 1 sequence", dest="a1", required=True)
parser.add_argument('-a2', help="adapter 2 sequence", dest="a2", required=True)
parser.add_argument('-maxns', help="max number of N's", dest="maxns", required=True)
parser.add_argument('-minq', help="minimum quality", dest="minq", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-s', help="stats file", dest="stats", required=True)
args = parser.parse_args()

read1i=args.read1i
read2i=args.read2i
read1o=args.read1o
read2o=args.read2o
a1=args.a1
a2=args.a2
maxns=args.maxns
minq=args.minq
threads=args.threads
stats=args.stats



# Run

statsfile=open(str(stats),"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
statsfile.write("Statistic\tValue \r\n".format(current_time))

#Get initial stats
reads = 0
bases = 0
#If gzipped
import os
if str(read1i).endswith('.gz'):
    with gzip.open(str(read1i), 'rb') as read:
        for id in read:
            seq = next(read)
            reads += 1
            bases += len(seq.strip())*2
            next(read)
            next(read)
else:
    with open(str(read1i), 'rb') as read:
        for id in read:
            try:
                seq = next(read)
                reads += 1
                bases += len(seq.strip())*2
                next(read)
                next(read)
            except:
                break
statsfile.write("Input reads\t{0} ({1} bases)\r\n".format(reads,bases))
statsfile.close()



# Run AdapterRemoval
if not os.path.exists(str(read1o)):
    if (a1 and a2):
        qualfiltCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.1.3 && AdapterRemoval --file1 '+read1i+' --file2 '+read2i+' --output1 '+read1o+' --output2 '+read2o+' --trimqualities --trimns --maxns '+maxns+' --minquality '+minq+' --threads '+threads+' --adapter1 '+a1+' --adapter2 '+a2+''
        subprocess.check_call(qualfiltCmd, shell=True)

    else: # default Illumina adapters will be used
        qualfiltCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.1.3 && AdapterRemoval --file1 '+read1i+' --file2 '+read2i+' --output1 '+read1o+' --output2 '+read2o+' --trimqualities --trimns --maxns '+maxns+' --minquality '+minq+' --threads '+threads+''
        subprocess.check_call(qualfiltCmd, shell=True)



#Get stats after quality filtering
reads = 0
bases = 0
with open(str(read1o), 'rb') as read:
    for id in read:
        try:
            seq = next(read)
            reads += 1
            bases += len(seq.strip())
            next(read)
            next(read)
        except:
            break

#Print stats to stats file
statsfile=open(str(str(stats)),"a+")
statsfile.write("Quality filtered reads\t{0} ({1} bases)\r\n".format(reads,bases))
statsfile.close()
