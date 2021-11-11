 #13.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import re
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-fq_path', help="path to .fastq files", dest="fq_path", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-obam_b', help="output bam file base", dest="obam_b", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-checkpoint', help="checkpoint file", dest="checkpoint", required=True)
args = parser.parse_args()


a=args.a
fq_path=args.fq_path
t=args.t
obam_b=args.obam_b
ID=args.ID
log=args.log
checkpoint=args.checkpoint


# Run

if os.path.exists(obam_b):
    pass

if not os.path.exists(obam_b):
    os.makedirs(obam_b)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as log:
        log.write('\t\t'+current_time+'\tCoAssembly Mapping step - '+ID+'\n')
        log.write('The original metagenomic reads are being mapped to the indexed assembly so coverage info can be retrieved.\n\n')


    # Get read1 and read2 paths
    reads1=glob.glob(fq_path+'/*_1.fastq*')

    for read1 in reads1:
        sampleID=os.path.basename(read1) # get sample ID from read2
        if sampleID.endswith('.gz'):
            sampleID=sampleID.replace('_1.fastq.gz','')
            read2=fq_path+'/'+sampleID+'_2.fastq.gz'
        else:
            sampleID=sampleID.replace('_1.fastq','')
            read2=fq_path+'/'+sampleID+'_2.fastq'

        obam=obam_b+'/'+sampleID+'.mapped.bam' # output bam path
        unmapped_r1=obam_b+'/'+sampleID+'_unmapped_1.fastq.gz'
        unmapped_r2=obam_b+'/'+sampleID+'_unmapped_2.fastq.gz'

        if not os.path.exists(str(obam)): # run bwa if output bam does not exist
            mappingCmd='module load tools samtools/1.11 bowtie2/2.4.2 \
            && bowtie2 \
            --time \
            --threads 10 \
            --rg-id "'+ID+'" \
            -x '+a+' \
            -1 '+read1+' \
            -2 '+read2+' \
            | samtools view -b -@ 10 - | samtools sort -@ 10 -T '+obam+'.'+sampleID+' -o '+obam+''
            # -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+a+' '+read1+' '+read2+' \
            # | samtools view -b - | samtools sort -T '+obam+'.'+sampleID+' -o '+obam+''
            subprocess.Popen(mappingCmd, shell=True).wait()

#        open('+checkpoint+', 'a').close()
        # if not os.path.exists(str(unmapped_r1)):
        #     unmappedCmd='module load tools samtools/1.11 && \
        #     samtools view -@ 10 \
        #     -b -f12 '+obam+' | \
        #     samtools fastq -@ 10 \
        #     -c 6 \
        #     -1 '+unmapped_r1+' \
        #     -2 '+unmapped_r1+' -'
        #     subprocess.Popen(unmappedCmd, shell=True).wait()
