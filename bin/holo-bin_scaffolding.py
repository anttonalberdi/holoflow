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


bin_dir= (bin_dir+'/dereplicated_genomes')
# Run
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tBin Scaffolding step - '+ID+'\n')
        logi.write('Scaffolds are build from the contigs found in every metagenomic bin by SSPACE.\n\n')


    binlist = glob.glob(str(bin_dir)+"/*.fa")
    for bin in binlist:
        bin_name=os.path.basename(bin)
        bin_name = bin_name.replace(".contigs.fa","")
        lib_file=str(out_dir+'/'+bin_name+'.lib')

        #Create library file
            # Insertion size between paired reads: 150
            # Maximum allowed error: 1
        libCmd='printf "'+ID+' bwa '+fq_dir+'/'+bin_name+'_1.fastq '+fq_dir+'/'+bin_name+'_2.fastq 150 1 FR" >> '+lib_file+' && cat '+lib_file+''
        subprocess.check_call(libCmd, shell=True)
        #Run SSPACE
        sspaceCmd ='cd '+out_dir+' && module load tools perl/5.24.0 sspace-standard/3.0 parallel/20190522 && SSPACE_Standard_v3.0.pl -l '+lib_file+' -s '+bin+' -x 1 -T '+threads+' -o 5 -m 16 -k 2 -n 10 -b '+bin_name+''
        subprocess.check_call(sspaceCmd, shell=True)


    #Rearrange outputs
    for bin in binlist:
        bin_name=os.path.basename(bin)
        bin_name = bin_name.replace(".contigs.fa","")
        faoutpCmd='cp '+out_dir+'/'+bin_name+'/'+bin_name+'.final.scaffolds.fasta '+out_dir+'/../'+bin_name+'.fna'
        subprocess.check_call(faoutpCmd, shell=True)
        infoutCmd='cp '+out_dir+'/'+bin_name+'/'+bin_name+'.summaryfile.txt '+out_dir+'/../'+bin_name+'.info'
        subprocess.check_call(infoutCmd, shell=True)
        rmCmd='rm '+out_dir+'/* && mv *.info *.fna '+out_dir+''
        subprocess.check_call(rmCmd, shell=True)
