#01.10.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-genomes_dir', help="input bin directory", dest="genomes_dir", required=True)
parser.add_argument('-div', help="diversity in PhyloPhlAn", dest="diversity", required=True)
parser.add_argument('-pip', help="PhyloPhlAn pipeline to be used", dest="pip", required=True)
parser.add_argument('-ph_db', help="genomes data base to be used by PhyloPhlAn", dest="ph_db", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()



genomes_dir=args.genomes_dir
diversity=args.diversity
pip=args.pip
ph_db=args.ph_db
out_dir=args.out_dir
sample=args.sample
log=args.log
threads=args.threads


# Run
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tMAG Phylogenetic assignation step - Sample '+sample+'\n')
        logi.write('\n\n')

        #Run PhyloPhlAn
        if pip == 'concatenation':
            pp_configCmd ='module load tools anaconda3/4.4.0 phylophlan/3.0 && cd '+out_dir+' && phylophlan_write_default_configs.sh && phylophlan -i '+genomes_dir+' -d '+ph_db+' --diversity '+diversity+' -f '+out_dir+'/supermatrix_nt.cfg && rm supermatrix_aa.cfg supertree_nt.cfg supertree_aa.cfg'
            subprocess.check_call(pp_configCmd, shell=True)

        if pip == 'tree':
            pp_configCmd ='module load tools anaconda3/4.4.0 phylophlan/3.0 && cd '+out_dir+' && phylophlan_write_default_configs.sh && phylophlan -i '+genomes_dir+' -d '+ph_db+' --diversity '+diversity+' -f '+out_dir+'/supertree_nt.cfg && rm supermatrix_aa.cfg supermatrix_nt.cfg supertree_aa.cfg'
            subprocess.check_call(pp_configCmd, shell=True)


    with open(str(log),'a+') as logf:
        logf.write('\t\t'+current_time+'\tMetagenomics analysis with Holoflow are completed for sample '+sample+'\n')
