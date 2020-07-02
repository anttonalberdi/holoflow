#01.07.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-bam', help="assembly mapped bam", dest="bam", required=True)
parser.add_argument('-dastool_bd', help="dastool bin directory", dest="dt_bd", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="main_out_dir", required=True)
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)

args = parser.parse_args()


a=args.a
bam=args.bam
dt_bd=args.dt_bd
main_out_dir=args.main_out_dir
sample=args.sample
log=args.log
threads=args.threads



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

condaenvCmd='conda activate /home/projects/ku-cbd/data/envs/refinem-0.1.1' # if doesn't work, source
subprocess.check_call(condaenvCmd, shell=True)


    ### Refinement based on genome properties

scaffold_statsCmd='refinem scaffold_stats -c '+threads+' --genome_ext fa '+assembly_file+' '+dt_bd+' '+main_out_dir+' '+bam+'' #assembly mapping bam / INTERSECT assembly
subprocess.check_call(scaffold_statsCmd, shell=True)

outliersCmd='refinem outliers '+main_out_dir+'/scaffold_stats.tsv '+main_out_dir+''
subprocess.check_call(outliersCmd, shell=True)

filter_binsCmd='refinem filter_bins --genome_ext fa  '+dt_bd+' '+main_out_dir+'/outliers.tsv '+main_out_dir+'/1_genomeproperties/'
subprocess.check_call(filter_binsCmd, shell=True)



    ### Refinement based on taxonomy

callgenesCmd='refinem call_genes -c 40 --genome_ext fa '+dt_bd+' '+main_out_dir+'/2_taxonomy/genes'
subprocess.check_call(callgenesCmd, shell=True)

txnprofileCmd='refinem taxon_profile -c 40 --tmpdir '+main_out_dir+'/2_taxonomy/tmp '+main_out_dir+'/2_taxonomy/genes '+main_out_dir+'/scaffold_stats.tsv /home/projects/ku-cbd/people/antalb/databases/RefineM/gtdb_r89_protein_db.2019-09-27.faa.dmnd /home/projects/ku-cbd/people/antalb/databases/RefineM/gtdb_r89_taxonomy.2019-09-27.tsv '+main_out_dir+'/2_taxonomy/'
subprocess.check_call(txnprofileCmd, shell=True)

txnfilterCmd='refinem taxon_filter -c 40 '+main_out_dir+'/2_taxonomy/ '+main_out_dir+'/2_taxonomy/taxon_filter.tsv'
subprocess.check_call(txnfilterCmd, shell=True)


#Refinement based on 16S genes
# mkdir ${workdir}/bin_refinement/3_16s
# mkdir ${workdir}/bin_refinement/4_finalbins
ssuerrCmd='refinem ssu_erroneous -c 40 --genome_ext fa '+main_out_dir+'/2_taxonomy '+main_out_dir+'/2_taxonomy /home/projects/ku-cbd/people/antalb/databases/RefineM/gtdb_r80_ssu_db.2018-01-18.fna /home/projects/ku-cbd/people/antalb/databases/RefineM/gtdb_r80_taxonomy.2017-12-15.tsv ${workdir}/bin_refinement/3_16s'
subprocess.check_call(ssuerrCmd, shell=True)

ssfilterCmd='refinem filter_bins --genome_ext fa '+main_out_dir+'/2_taxonomy '+main_out_dir+'/3_16s/ssu_erroneous.tsv '+main_out_dir+'/4_finalbins && rm '+main_out_dir+'/4_finalbins/refinem.log'
subprocess.check_call(ssfilterCmd, shell=True)


with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMetagenomics analysis with Holoflow are completed for sample '+sample+'\n')
