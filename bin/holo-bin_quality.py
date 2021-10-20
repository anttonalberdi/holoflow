#02.11.2020

import subprocess
import argparse
import time
import sys
import os


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bin_dir', help="drep bin directory", dest="bin_dir", required=True)
parser.add_argument('-cov_file', help="coverage data file ", dest="cov_file", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()



bin_dir=args.bin_dir
cov_file=args.cov_file
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
input_checkM_table = bin_dir+'/Widb.csv'
if not os.path.isfile(out_dir+'/'+ID+'_binQuality.txt'):
    checkmCmd = 'module load anaconda2/4.0.0 hmmer/3.2.1 prodigal/2.6.3 pplacer/1.1.alpha17 \
    && checkm lineage_wf -t '+threads+' -x fa '+bin_dir+' '+out_dir+' -f '+out_dir+'/'+ID+'_binQuality.txt'
    subprocess.Popen(checkmCmd,shell=True).wait()

    rearrangeoutCmd =' sed -i "s/--//g" '+out_dir+'/'+ID+'_binQuality.txt \
    && sed -i "s/ \+ /\t/g" '+out_dir+'/'+ID+'_binQuality.txt'
    subprocess.Popen(rearrangeoutCmd,shell=True).wait()

# Plot quality - coverage
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)

if os.path.isfile(out_dir+'/'+ID+'_binQuality.txt'):
    plotCmd = 'module load tools gcc/5.4.0 intel/compiler/64/2018_update2 R/3.5.3-ICC-MKL \
    && Rscript '+curr_dir+'/holo-bin_quality.plot.R -cov_data '+str(cov_file)+' -qual_data '+out_dir+'/'+ID+'_binQuality.txt -ID '+ID+' -out_path '+out_dir+''
    subprocess.Popen(plotCmd,shell=True).wait()


# Run summary table
input_drep_table = bin_dir+'/final_bins_Info.csv'
input_checkM_table = bin_dir+'/Widb.csv'
summary_table_tmp = out_dir+'/'+ID+'_binQuality_Info.tmp.csv'
mag_table = out_dir+'/'+ID+'_binQuality_detail_Info.csv'
summary_table = out_dir+'/'+ID+'_binQuality_general_Info.csv'


if os.path.isfile(out_dir+'/'+ID+'_binQuality.txt'):
    summaryCmd = 'bash '+curr_dir+'/holo-bin_quality_table.sh '+input_drep_table+' '+input_checkM_table+' '+summary_table_tmp+' '+mag_table+' '+summary_table+''
    subprocess.Popen(summaryCmd,shell=True).wait()
