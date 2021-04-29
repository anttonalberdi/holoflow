#22.11.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time
import re
import numpy as np


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-fq_dir', help="input .fq directory", dest="fq_dir", required=True)
parser.add_argument('-bin_dir', help="input bin directory", dest="bin_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


fq_dir=args.fq_dir
bin_dir=args.bin_dir
out_dir=args.out_dir
ID=args.ID
log=args.log
threads=args.threads


# Run
# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\t\t'+current_time+'\tMAG Mapping step - '+ID+'\n')
    logi.write('MAGs are being mapped to the original metagenomic read files to assess its coverage.\n\n')


# Create MAGs file --> competitive mapping for each sample
mag_catalogue_file=out_dir+'/'+ID+'_MAG_Catalogue.fa'

if not (os.path.isfile(str(mag_catalogue_file))):
    with open(mag_catalogue_file,'w+') as magcat:

        maglist = glob.glob(str(bin_dir)+"/*.fa")
        for mag in maglist:
            mag_name=os.path.basename(mag)
            mag_name = mag_name.replace(".fa","")

            with open(mag,'r') as mag_data:
                for line in mag_data.readlines():
                    if line.startswith('>'):
                        line=line.replace('>','>'+mag_name+'-')
                        magcat.write(line)
                    else:
                        magcat.write(line)


# Index MAG catalogue file
IDXmag_catalogue_file=out_dir+'/'+ID+'_MAG_Catalogue.fa.fai'

if not (os.path.isfile(str(IDXmag_catalogue_file))):
    idxsamCmd='module load tools samtools/1.11 && samtools faidx '+mag_catalogue_file+''
    idxbwaCmd='module load tools bwa/0.7.15 && bwa index '+mag_catalogue_file+''

    #subprocess.Popen(idxbwaCmd, shell=True).wait()
    #subprocess.Popen(idxsamCmd, shell=True).wait()


# Initialize stats
stats_file = out_dir+'/'+ID+'.MAG_mapping_stats.txt'
sample_list = list()
mapped_reads_tmp = out_dir+'/'+ID+'.tmp_mapped.reads.txt'
total_reads_tmp = out_dir+'/'+ID+'.tmp_total.reads.txt'

if (os.path.isfile(str(IDXmag_catalogue_file))):
    readlist = glob.glob(str(fq_dir)+"/*.fastq*")
    samples = list()
    for file in readlist:
        read_name=''
        read_name=os.path.basename(file)
        if file.endswith('.gz'):
            extension = '.gz'
            read_name = re.sub('_[0-9]\.fastq.gz','',read_name)
        else:
            extension = ''
            read_name = re.sub('_[0-9]\.fastq','',read_name)
        samples.append(read_name)
    sample_list = sorted(set(samples))

    for sample in sample_list:
        # Map every sample to mag catalogue file (competitive mapping) - get one bam for every sample
        out_bam = out_dir+'/'+sample+'.bam'

        if extension == '.gz':
            read1 = fq_dir+'/'+sample+'_1.fastq.gz'
            read2 = fq_dir+'/'+sample+'_2.fastq.gz'
        else:
            read1 = fq_dir+'/'+sample+'_1.fastq'
            read2 = fq_dir+'/'+sample+'_2.fastq'

        mapbinCmd='module load tools samtools/1.11 bwa/0.7.15 && bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+mag_catalogue_file+' '+read1+' '+read2+' | samtools view -b - | samtools sort -T '+out_dir+'/'+ID+' -o '+out_bam+''
        #subprocess.Popen(mapbinCmd, shell=True).wait()

        # extract not-mapped to the reference genome reads + keep reference bam
        read1_not=out_dir+'/'+sample+'_notMAGmap_1.fastq.gz'
        read2_not=out_dir+'/'+sample+'_notMAGmap_2.fastq.gz'
        refbamCmd = 'module load tools samtools/1.11 && samtools view -T '+mag_catalogue_file+' -b -f12 '+out_bam+' | samtools fastq -1 '+read1_not+' -2 '+read2_not+' -'
        subprocess.Popen(refbamCmd, shell=True).wait()


######################## Stats ########################

        # Get total number of initial reads bases
        # samtools view -c
        totalCmd='module load tools samtools/1.11 && samtools view -c '+out_bam+' >> '+total_reads_tmp+''
        #subprocess.Popen(totalCmd, shell=True).wait()


        # Get mapped number of reads
        # samtools view -c -F 4
        mappedCmd='module load tools samtools/1.11 && samtools view -c -F 4 '+out_bam+' >> '+mapped_reads_tmp+''
        #subprocess.Popen(mappedCmd, shell=True).wait()


    ## Build stats file
    # Write sample IDs
    stats = open(stats_file,'w+')
    sample_list.insert(0,'Sample_ID')
    stats.write(('\t').join(sample_list)+'\n')

        # Retrieve all numbers of MAPPED reads
    with open(mapped_reads_tmp,'r+') as mapped_reads_file:
        mapped_reads = list()
        for line in mapped_reads_file.readlines():
            mapped_reads.append(line.strip())
    #os.remove(mapped_reads_tmp)

        # Retrieve all numbers of TOTAL reads
    with open(total_reads_tmp,'r+') as total_reads_file:
        total_reads = list()
        for line in total_reads_file.readlines():
            total_reads.append(line.strip())
    #os.remove(total_reads_tmp)


    # Write number of mapped reads per sample
    #stats.write('Mapped Reads'+'\t'+('\t').join(mapped_reads)+'\n')

        # Calculate percentage of mapped reads from: (mapped reads/ total reads) * 100
    mapped_reads = np.array(mapped_reads).astype(int)
    total_reads = np.array(total_reads).astype(int)
    percentages = np.divide(mapped_reads,total_reads)
    percentages = (percentages*100)
    percentages = percentages.round(decimals=2).tolist() # true division

    # Write percentagesfinal_tips = (',').join('"{0}"'.format(tip) for tip in final_tips)
    #stats.write('% Mapped Reads'+'\t'+('\t').join(str(perc) for perc in percentages))
