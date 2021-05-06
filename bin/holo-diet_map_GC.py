#06.05.2021 - Holoflow 0.1.
import subprocess
import argparse
import os
import time
import glob


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-fna', help="nucleotidic sequences predicted ORFs", dest="fna", required=True)
parser.add_argument('-fq_dir', help="unmapped reads to MAGs fq directory", dest="fq_dir", required=True)
parser.add_argument('-out_dir', help="out_dir", dest="out_dir", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


fna=args.fna
fq_dir=args.fq_dir
out_dir=args.out_dir
t=args.threads
ID=args.ID
log=args.log



# Run
# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\tHOLOFLOW\tMETAGENOMICS\n\t\t'+current_time+'\t - '+ID+'\n')
    logi.write('   \n\n')

# index gene catalogue file
if not os.path.exists(fna+'.fai'):
    idxsamCmd='module load tools samtools/1.11 && samtools faidx '+fna+''
    idxbwaCmd='module load tools bwa/0.7.15 && bwa index '+fna+''

    subprocess.Popen(idxbwaCmd, shell=True).wait()
    subprocess.Popen(idxsamCmd, shell=True).wait()



if os.path.exists(fna+'.amb'):
# Get read1 and read2 paths
    reads1=glob.glob(fq_dir+'/*_1.fastq.gz')

    for read1 in reads1:
        sampleID=os.path.basename(read1)
        sampleID=sampleID.replace('_1.fastq.gz','')

        read2=fq_dir+'/'+sampleID+'_2.fastq.gz'
        obam=obam_b+'/'+ID+'.'+sampleID+'.MAG_unmapped.bam'

        if not os.path.exists(str(obam)):
            mappingCmd='module load tools samtools/1.11 bwa/0.7.15 && bwa mem -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+fna+' '+read1+' '+read2+' | samtools view -b - | samtools sort -T '+obam+'.'+sampleID+' -o '+obam+''
            subprocess.Popen(mappingCmd, shell=True).wait()
