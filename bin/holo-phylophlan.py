#01.10.2020 - Holoflow 0.1.
################################### NOT IN USE NOW ##################################

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
parser.add_argument('-ssp', help="SSPACE used or not", dest="ssp", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()



genomes_dir=args.genomes_dir
diversity=args.diversity
pip=args.pip
ph_db=args.ph_db
out_dir=args.out_dir
ssp=args.ssp
ID=args.ID
log=args.log
threads=args.threads


# Run
# if not (os.path.exists(str(out_dir))):
#     os.mkdir(str(out_dir))

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\t\t'+current_time+'\tMAG Phylogenetic assignation step - '+ID+'\n')
    logi.write('\n\n')

if not (ssp): #drep output files have .fa extension, PhyloPhlAn requires .fna for nucl.
    genomelist=glob.glob(str(genomes_dir)+"/*.fa")
    for genome in genomelist:
        genome_n=genome.replace(".fa",".fna")
        genomeCmd='mv '+genome+' '+genome_n+''
        subprocess.check_call(genomeCmd,shell=True)


#Run PhyloPhlAn
if pip == 'concatenation':
    pp_configCmd ='module load tools anaconda3/4.4.0 phylophlan/3.0 && cd '+out_dir+'/.. && phylophlan_write_config_file -o holoflow_matrix_config_nt.cfg -d a --force_nucleotides --db_aa diamond --map_aa diamond --map_dna diamond --msa muscle --tree1 fasttree && phylophlan -i '+genomes_dir+' -d '+ph_db+' --diversity '+diversity+' --force_nucleotides -f '+out_dir+'/../holoflow_matrix_config_nt.cfg -o Matrix_Database'
    subprocess.check_call(pp_configCmd, shell=True)

if pip == 'tree':
    pp_configCmd ='module load tools anaconda3/4.4.0 phylophlan/3.0 && cd '+out_dir+'/.. && phylophlan_write_config_file -o holoflow_tree_config_nt.cfg -d a --force_nucleotides --db_aa diamond --map_aa diamond --map_dna diamond --msa muscle --tree1 fasttree --gene_tree1 fasttree --gene_tree2 ramxl &&  phylophlan -i '+genomes_dir+' -d '+ph_db+' --diversity '+diversity+' --force_nucleotides -f '+out_dir+'/../holoflow_tree_config_nt.cfg -o Tree_Database --maas /services/tools/phylophlan/3.0/lib/python3.7/site-packages/phylophlan/phylophlan_substitution_models/phylophlan.tsv'
    subprocess.check_call(pp_configCmd, shell=True)


with open(str(log),'a+') as logf:
    logf.write('\t\t'+current_time+'\tMetagenomics analysis with Holoflow are completed for ID '+ID+'\n')
