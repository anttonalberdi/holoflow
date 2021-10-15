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
    logi.write('The reads not included in the MAG set are mapped to the gene catalogue created by Prodigal 2.6.3.\n\n')


# index gene catalogue file: .fna predicted sequences by prodigal
if not os.path.exists(fna+'.fai'):
    idxsamCmd='module load tools samtools/1.11 && samtools faidx '+fna+''
    idxbt2Cmd='module load tools bowtie2/2.4.2 && bowtie2-build --threads '+t+' '+fna+' '+fna+''

    subprocess.Popen(idxbt2Cmd, shell=True).wait()
    subprocess.Popen(idxsamCmd, shell=True).wait()



if os.path.exists(fna+'.rev.2.bt2l'):
# Get read1 and read2 paths #### reads that were not mapped to MAGs
    reads1=glob.glob(fq_dir+'/*_1.fastq*')


    for read1 in reads1:
        sampleID=os.path.basename(read1)
        sampleID=sampleID.replace('_1.fastq.gz','')

        read2=fq_dir+'/'+sampleID+'_2.fastq.gz'
        obam=out_dir+'/'+ID+'.'+sampleID+'.MAG_unmapped.bam' # output bam

        if not os.path.exists(out_dir):
            mkdirCmd='mkdir -p '+out_dir+''
            subprocess.Popen(mkdirCmd,shell=True).wait()

        if not os.path.exists(str(obam)):   # run mapping
            mappingCmd='module load tools samtools/1.11 bowtie2/2.4.2 && \
            bowtie2 \
            --threads '+t+' \
            --rg-id "'+ID+'" \
            -x '+fna+' \
            -1 '+read1+' \
            -2 '+read2+' \
            | samtools view -@ '+t+' -b - | samtools sort -@ '+t+' -T '+obam+'.'+sampleID+' -o '+obam+''
            subprocess.Popen(mappingCmd, shell=True).wait()
