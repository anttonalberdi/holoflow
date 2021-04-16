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
parser.add_argument('-msep', help="mate separator between 1,2 reads", dest="msep", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
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
msep=args.msep
log=args.log
threads=args.threads
stats=args.stats



# Run
statsfile=open(str(stats),"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
statsfile.write("Statistic\tValue \r\n".format(current_time))

if (os.path.exists(read1i)):
    compressCmd1='gunzip '+read1i+' '+read2i+''
    subprocess.Popen(compressCmd1,shell=True).wait()
    read1i = read1i.replace('.gz','')
    read2i = read2i.replace('.gz','')


#Get initial stats
reads = 0
bases = 0
#If gzipped
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


# Write to log
with open(str(log),'a+') as log:
    log.write('\tHOLOFLOW\tPREPROCESSING\n\t\t'+current_time+'\tQuality Filtering step\n')
    log.write('Those reads with a minimum quality of '+minq+' are being removed.\nThe sequencing adapters of all reads as well.\n\n')




# Run AdapterRemoval
# output --gzip files
# use a diferent separator of reads 
if not (msep == "default"):
    if not os.path.exists(str(read1o)):
        # different adapters than default
        if not ((a1 == "default") and (a2 == "default")):
            qualfiltCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.2.4 && AdapterRemoval --file1 '+read1i+' --file2 '+read2i+' --mate-separator '+msep+' --output1 '+read1o+' --output2 '+read2o+' --gzip --trimqualities --trimns --maxns '+maxns+' --minquality '+minq+' --threads '+threads+' --adapter1 '+a1+' --adapter2 '+a2+''
            subprocess.check_call(qualfiltCmd, shell=True)

        else: # default Illumina adapters will be used
            qualfiltCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.2.4 && AdapterRemoval --file1 '+read1i+' --file2 '+read2i+' --mate-separator '+msep+' --output1 '+read1o+' --output2 '+read2o+' --gzip --trimqualities --trimns --maxns '+maxns+' --minquality '+minq+' --threads '+threads+''
            subprocess.check_call(qualfiltCmd, shell=True)
else:
    if not os.path.exists(str(read1o)):
        if not ((a1 == "default") and (a2 == "default")):
            qualfiltCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.2.4 && AdapterRemoval --file1 '+read1i+' --file2 '+read2i+' --output1 '+read1o+' --output2 '+read2o+' --gzip --trimqualities --trimns --maxns '+maxns+' --minquality '+minq+' --threads '+threads+' --adapter1 '+a1+' --adapter2 '+a2+''
            subprocess.check_call(qualfiltCmd, shell=True)

        else: # default Illumina adapters will be used
            qualfiltCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.2.4 && AdapterRemoval --file1 '+read1i+' --file2 '+read2i+' --output1 '+read1o+' --output2 '+read2o+' --gzip --trimqualities --trimns --maxns '+maxns+' --minquality '+minq+' --threads '+threads+''
            subprocess.check_call(qualfiltCmd, shell=True)



#Get stats after quality filtering
# read --gzip files
reads = 0
bases = 0
with gzip.open(str(read1o), 'rt') as read:
    for id in read:
        try:
            seq = next(read)
            reads += 1
            bases += len(seq.strip())*2
            next(read)
            next(read)
        except:
            break

# re-compress inputs
if (os.path.exists(read1o)):
    compressCmd2='gzip '+read1i+' '+read2i+''
    subprocess.Popen(compressCmd2,shell=True).wait()

#Print stats to stats file
statsfile=open(str(str(stats)),"a+")
statsfile.write("Quality filtered reads\t{0} ({1} bases)\r\n".format(reads,bases))
statsfile.close()
