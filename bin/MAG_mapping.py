#24.09.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


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
        logi.write('\t\t'+current_time+'\tBin Mapping step - '+ID+'\n')
        logi.write('MAGs are being mapped to the original metagenomic read files to assess its coverage.\n\n')


    binlist = glob.glob(str(bin_dir)+"/*.fa")
    for bin in binlist:
        bin_name=os.path.basename(bin)
        bin_name = bin_name.replace(".contigs.fa","")
        lib_file=str(out_dir+'/'+bin_name+'.lib')

        #Create library file
            # Insertion size between paired reads: 150
            # Maximum allowed error: 1
        libCmd='module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+t+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:ID" '+bin+' '+fq_dir+'/'+bin_name+'_1.fastq '+fq_dir+'/'+bin_name+'_2.fastq | samtools view -b - | samtools sort - > '+obam+''
        subprocess.check_call(libCmd, shell=True)
