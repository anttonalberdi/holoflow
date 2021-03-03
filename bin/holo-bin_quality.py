#02.11.2020

import subprocess
import argparse
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bin_dir', help="drep bin directory", dest="bin_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()



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
        logi.write('\t\t'+current_time+'\tBin Quality step - '+ID+'\n')
        logi.write('\n\n')


    ## RUN

    bin_dir=bin_dir+'/dereplicated_genomes'

    checkmCmd = 'module load anaconda2/4.0.0 hmmer/3.2.1  prodigal/2.6.3 pplacer/1.1.alpha17 && checkm lineage_wf -t '+threads+' -x fa '+bin_dir+' '+out_dir+' -f '+out_dir+'/'+ID+'_binQuality.txt'
    subprocess.Popen(checkmCmd,shell=True).wait()
