## 26.02.21 - Holoflow 0.1
import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-filt_dir', help="filtered variants directory", dest="filt_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-gmap', help="gmap", dest="gmap", required=True)
args = parser.parse_args()


filt_dir=args.filt_dir
out_dir=args.out_dir
chr_list=args.chr_list
ID=args.ID
log=args.log
threads=args.threads
gmap=args.gmap


## Run
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tPhasing of HD data - '+ID+'\n')
        logi.write(' \n\n')

    chromosome_list = list()
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            chromosome_list.append(chr.strip())

    for CHR in chromosome_list:
        input = filt_dir+'/'+ID+'.HD_SNPs_'+CHR+'.vcf.gz'
        output = out_dir+'/'+ID+'_'+CHR+'.filt_phased.vcf.gz'


        if not (gmap == 'False'):
            phasingCmd= 'module load shapeit4/4.1.3 && shapeit4 --input '+input+' --map '+gmap+' --region '+CHR+' --thread '+threads+' --output '+output+' --sequencing'
            subprocess.Popen(phasingCmd,shell=True).wait()

        else:
            phasingCmd= 'module load shapeit4/4.1.3 && shapeit4 --input '+input+' --region '+CHR+' --thread '+threads+' --output '+output+' --sequencing'
            subprocess.Popen(phasingCmd,shell=True).wait()

        # Index phased panel
        idxCmd='module load tabix/1.2.1 && tabix '+output+''
        subprocess.Popen(idxCmd,shell=True).wait()
