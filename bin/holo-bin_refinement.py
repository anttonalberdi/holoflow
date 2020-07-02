#01.07.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bin_dir', help="assembly file", dest="a", required=True)
parser.add_argument('-out_dir', help="metabat bin table", dest="bt_mtb", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

bin_dir=args.bin_dir
out_dir=args.out_dir
sample=args.sample
log=args.log



# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tRefineM Bin Refinement step - Sample '+sample+'\n')
    log.write('Based on genome properties and taxonomy, RefineM will take the Dastool bins merged from Maxbin and Metabat2\nand try to increase its completeness while reducing the redundancy. \n\n')



# Filter assembly file - only those contigs in dastool




# RefineM
refinemDependenciesCmd='module load tools anaconda3/4.4.0 kronatools/2.7 diamond/0.9.29'
subprocess.check_call(refinemDependenciesCmd, shell=True)

source activate /home/projects/ku-cbd/data/envs/refinem-0.1.1 # activate conda environment - HOW here?


    ### Refinement based on genome properties

scaffold_statsCmd='refinem scaffold_stats -c '+threads+' --genome_ext fa '+assembly_file+' '+dastool_bins_dir+' '+main_output_dir+' '+bam_output+''
subprocess.check_call(scaffold_statsCmd, shell=True)

outliersCmd='refinem outliers '+main_output_dir+'/scaffold_stats.tsv '+main_output_dir+''
subprocess.check_call(outliersCmd, shell=True)

filter_binsCmd='refinem filter_bins --genome_ext fa  '+dastool_bins_dir+' '+main_output_dir+'/outliers.tsv '+main_output_dir+'/1_genomeproperties/'
subprocess.check_call(filter_binsCmd, shell=True)



    ### Refinement based on taxonomy

callgenesCmd='refinem call_genes -c 40 --genome_ext fa '+dastool_bins_dir+' '+main_output_dir+'/2_taxonomy/genes'
subprocess.check_call(callgenesCmd, shell=True)

txnprofileCmd='refinem taxon_profile -c 40 --tmpdir '+main_output_dir+'/2_taxonomy/tmp '+main_output_dir+'/2_taxonomy/genes '+main_output_dir+'/scaffold_stats.tsv /home/projects/ku-cbd/people/antalb/databases/RefineM/gtdb_r89_protein_db.2019-09-27.faa.dmnd /home/projects/ku-cbd/people/antalb/databases/RefineM/gtdb_r89_taxonomy.2019-09-27.tsv '+main_output_dir+'/2_taxonomy/'
subprocess.check_call(txnprofileCmd, shell=True)

txnfilterCmd='refinem taxon_filter -c 40 '+main_output_dir+'/2_taxonomy/ '+main_output_dir+'/2_taxonomy/taxon_filter.tsv'
subprocess.check_call(txnfilterCmd, shell=True)



with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMetagenomics analysis with Holoflow are completed for sample '+sample+'\n')
