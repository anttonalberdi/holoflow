#27.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-bt_mtb', help="metabat bin table", dest="bt_mtb", required=True)
parser.add_argument('-bt_mxb', help="maxbin bin table", dest="bt_mxb", required=True)
parser.add_argument('-p', help="prodigal predicted proteins", dest="p", required=True)
parser.add_argument('-o', help="output main dir", dest="o", required=True)
parser.add_argument('-se', help="search engine", dest="se", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-db', help="dastool database directory", dest="db", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

a=args.a
bt_mtb=args.bt_mtb
bt_mxb=args.bt_mxb
p=args.p
o=args.o
se=args.se
t=args.t
db=args.db
sample=args.sample
log=args.log



# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tDASTool Bin Refinement step - Sample '+sample+'\n')
    log.write('The binning results from MaxBin and Metabat2 are integrated by DASTool to produce one only non-redundant\nset of bins between them.\n\n')


dastoolDependencies='module unload maxbin/2.2.7 fraggenescan/1.31 perl/5.20.2 && module load tools gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ prodigal/2.6.3 das_tool/1.1.1 diamond/0.9.24 usearch/11.0.667'
dastoolCmd=''+dastoolDependencies+' && DAS_Tool -i '+bt_mxb+','+bt_mtb+' -c '+a+' -o '+o+' --proteins '+p+' -l maxbin,metabat --search_engine '+se+' -t '+t+' --db_directory '+db+' --write_bins 1'
subprocess.check_call(dastoolCmd, shell=True)


# Move definitive bins to final directory
binfiles = glob.glob(os.path.join(str(o),'*.fa'))
for b in binfiles:
    shutil.move(b, str(''+o+'.bin'))


# Add relevant info to log
with open(str(log),'a+') as log:
    log.write('\t\tDASTool MaxBin bins evaluation - Sample '+sample+'\n')
    with open(str(''+o+'_maxbin.eval'),'r') as mxb_eval:
        log.write(''+mxb_eval+'\n')
    log.write('\t\tDASTool Metabat2 bins evaluation - Sample '+sample+'\n')
    with open(str(''+o+'_metabat.eval'),'r') as mtb_eval:
        log.write(''+mtb_eval+'\n')
    log.write('\t\tDASTool Bin Merging Summary - Sample '+sample+'\n')
    with open(str(''+o+'_DASTool_summary.txt'),'r') as summary:
        log.write(''+summary+'\n\n')


mvinfoCmd='mv '+o+'_maxbin.eval '+o+'_metabat.eval '+o+'_DASTool_summary.txt ..'
subprocess.check_call(mvinfoCmd, shell=True)


with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMetagenomics analysis with Holoflow are completed for sample '+sample+'\n')
