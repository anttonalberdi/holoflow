#22.11.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-mag_dir', help="input mapped MAGs to .fastq directory", dest="mag_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


mag_dir=args.mag_dir
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
        logi.write('\t\t'+current_time+'\tMAG Coverage step - '+ID+'\n')
        logi.write('\n\n')

    # Extract MAGs coverage from bam files
    mapped_list = glob.glob(str(mag_dir)+"/*.bam")
    for bam in mapped_list:
        sample=''
        sample=os.path.basename(bam)
        sample=sample.replace(".bam","")
        covCmd='module load tools bedtools/2.28.0 && bedtools genomecov -ibam '+bam+' > '+out_dir+'/'+sample+'_MAGcoverage.bed'
        subprocess.Popen(covCmd, shell=True).wait()
