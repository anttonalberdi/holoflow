#24.09.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time
import re


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
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

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
                mag_name = mag_name.replace(".contigs.fa","")

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
        idxsamCmd='module load tools samtools/1.9 && samtools faidx '+mag_catalogue_file+''
        idxbwaCmd='module load tools bwa/0.7.15 && bwa index '+mag_catalogue_file+''

        subprocess.Popen(idxbwaCmd, shell=True).wait()
        subprocess.Popen(idxsamCmd, shell=True).wait()


    if (os.path.isfile(str(IDXmag_catalogue_file))):
        readlist = glob.glob(str(fq_dir)+"/*.fastq")
        samples = list()
        for file in readlist:
            read_name=''
            read_name=os.path.basename(file)
            read_name = re.sub('_[0-9]\.fastq','',read_name)
            samples.append(read_name)

        sample_list = set(samples)
        for sample in sample_list:
            # Map every sample to mag catalogue file (competitive mapping) - get one bam for every sample
            out_bam=out_dir+'/'+sample+'.bam'
            mapbinCmd='module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+mag_catalogue_file+' '+fq_dir+'/'+sample+'_1.fastq '+fq_dir+'/'+sample+'_2.fastq | samtools view -b - | samtools sort - > '+out_bam+''
            subprocess.Popen(mapbinCmd, shell=True).wait()
