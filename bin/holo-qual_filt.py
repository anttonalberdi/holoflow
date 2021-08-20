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
parser.add_argument('-maxns', help="max number of N's, default 5", dest="maxns", required=True)
parser.add_argument('-minq', help="minimum quality PHRED, default 30 -> 0.1%", dest="minq", required=True)
parser.add_argument('-minlen', help="minimum length of reads, default 35", dest="minlen", required=True)
parser.add_argument('-lowcomplexfilt', help="Enable low complexity read filtering", dest="-lowcomplexfilt", required=True)
parser.add_argument('-complexthreshold', help="Threhold for complexity, default 30%", dest="-complexthreshold", required=True)
#parser.add_argument('-msep', help="mate separator between 1,2 reads", dest="msep", required=True)
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
minlen=args.minlen
lowcomplexfilt=args.lowcomplexfilt
complexthreshold=args.complexthreshold
#msep=args.msep
log=args.log
threads=args.threads
stats=args.stats



# Run

# write to stats
statsfile=open(str(stats),"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
statsfile.write("Statistic\tValue \r\n".format(current_time))

### Raph: Removing this step, as it takes time, and we can get this information
### from the output of fastp -- which is quicker than decompressing.
### Keeping .stats file generation, as it's currently needed in the snakemake file.

#Get initial stats
#reads = 0
#bases = 0
#If gzipped
#with gzip.open(str(read1i), 'rt') as read:
#    for id in read:
#        try:
#            seq = next(read)
#            reads += 1
#            bases += len(seq.strip())*2
#            next(read)
#            next(read)
#        except:
#            break
#statsfile.write("Input reads\t{0} ({1} bases)\r\n".format(reads,bases))
#statsfile.close()


# Write to log
with open(str(log),'a+') as log:
    log.write('\tHOLOFLOW\tPREPROCESSING\n\t\t'+current_time+'\tQuality Filtering step\n')
    log.write('Reads with a minimum quality of '+minq+' are being removed.\n
               Trimming sequencing adapters of all reads as well.\n\n')



# Run Fastp
### Raph: Removing custom mate separator option (seems niche?)
# use a diferent separator of reads
if not (lowcomplexfilt == "true"):
    if not os.path.exists(str(read1o)):
        # different adapters than default
        if not ((a1 == "default") and (a2 == "default")):
            qualfiltCmd = 'module unload gcc tools ngs \
            && module load tools fastp/0.20.1 \
            && fastp \
            --in1 '+read1i+' --in2 '+read2i+' \
            --out1 '+read1o+' --out2 '+read2o+' \
            --compression \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit '+maxns+' \
            --qualified_quality_phred '+minq+' \
            --length_required '+minlen+'\
            --overrepresentation_analysis \
            --thread '+threads+' \
            --adapter_sequence '+a1+' \
            --adapter_sequence_r2 '+a2+''
            subprocess.check_call(qualfiltCmd, shell=True)

        else: # default Illumina adapters will be used
            qualfiltCmd = 'module unload gcc tools ngs \
            && module load tools fastp/0.20.1 \
            && fastp \
            --in1 '+read1i+' --in2 '+read2i+' \
            --out1 '+read1o+' --out2 '+read2o+' \
            --compression \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit '+maxns+' \
            --qualified_quality_phred '+minq+' \
            --length_required '+minlen+'\
            --overrepresentation_analysis \
            --thread '+threads+' \
            subprocess.check_call(qualfiltCmd, shell=True)
else:
    if not os.path.exists(str(read1o)):
        if not ((a1 == "default") and (a2 == "default")):
            qualfiltCmd = 'module unload gcc tools ngs \
            && module load tools fastp/0.20.1 \
            && fastp \
            --in1 '+read1i+' --in2 '+read2i+' \
            --out1 '+read1o+' --out2 '+read2o+' \
            --compression \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit '+maxns+' \
            --qualified_quality_phred '+minq+' \
            --length_required '+minlen+'\
            --low_complexity_filter \
            --complexity_threshold '+complexthreshold+'\
            --overrepresentation_analysis \
            --thread '+threads+' \
            --adapter_sequence '+a1+' \
            --adapter_sequence_r2 '+a2+''
            subprocess.check_call(qualfiltCmd, shell=True)

        else: # default Illumina adapters will be used
            qualfiltCmd = 'module unload gcc tools ngs \
            && module load tools fastp/0.20.1 \
            && fastp \
            --in1 '+read1i+' --in2 '+read2i+' \
            --out1 '+read1o+' --out2 '+read2o+' \
            --compression \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit '+maxns+' \
            --qualified_quality_phred '+minq+' \
            --length_required '+minlen+'\
            --low_complexity_filter \
            --complexity_threshold '+complexthreshold+'\
            --overrepresentation_analysis \
            --thread '+threads+' \
            subprocess.check_call(qualfiltCmd, shell=True)



#Get stats after quality filtering
# read --gzip files
#reads = 0
#bases = 0
#with gzip.open(str(read1o), 'rt') as read:
#    for id in read:
#        try:
#            seq = next(read)
#            reads += 1
#            bases += len(seq.strip())*2
#            next(read)
#            next(read)
#        except:
#            break


#Print stats to stats file
statsfile=open(str(str(stats)),"a+")
#statsfile.write("Quality filtered reads\t{0} ({1} bases)\r\n".format(reads,bases))
statsfile.close()
