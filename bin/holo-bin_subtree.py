#19.01.21
import subprocess
import argparse
import os
import sys
import glob
import time
import re


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-tree_dir', help="gtdbtk phylogenetic trees directory", dest="tree_dir", required=True)
parser.add_argument('-bin_dir', help="dereplicated bins dir", dest="bin_dir", required=True)
parser.add_argument('-bac_o', help="output BAC subtree", dest="bac_o", required=True)
parser.add_argument('-ar_o', help="output AR subtree", dest="ar_o", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


tree_dir=args.tree_dir
bin_dir=args.bin_dir
bac_o=args.bac_o
ar_o=args.ar_o
ID=args.ID
log=args.log


# Run
#if not (os.path.isfile(bac_o)):

# Define in and out tree paths
in_paths = sorted(glob.glob(tree_dir+'/*.tree'))
out_paths = [ar_o,bac_o]

# In case bins come from individually assembled samples: get all sample IDs in group
# If bins come from coassembly, only one ID will be in the list
bin_names = [os.path.basename(x) for x in glob.glob(bin_dir+'/*.fa')]
ID_list = list()
for ID in bin_names:
    ID = ID.split('.')[0]   # Assume user won't use '.' in sample ID :()
    if not ID in ID_list:
        ID_list.append(ID)


##### Subtract group's tree tips - omit gtdbtk's entries
for i in range(len(in_paths)):
    tree_path = in_paths[i]
    out_tree_path = out_paths[i]
    tree_data = str()
    sample_tips = list()

    # Read in tree
    with open(tree_path,'r+') as tree:
        for line in tree.readlines():
            tree_data+=line

    # Pattern search for user's bins
    for ID in ID_list:
        # Find between 1 and unlimited case insensitive letters (ID), this can include numbers or not.
        # After that a . followed by three lower-case letters (mtb,cct,mxb) followed by 1,2,3 or 4 numbers (binner bin number)
        # followed by ".fa"
        match = re.findall(str(ID)+'[a-zA-Z]+[0-9]*\.{1}[a-z]{3}[0-9]{1}|[a-zA-Z]+[0-9]*\.{1}[a-z]{3}[0-9]{2}|[a-zA-Z]+[0-9]*\.{1}[a-z]{3}[0-9]{3}|[a-zA-Z]+[0-9]*\.{1}[a-z]{3}[0-9]{4}',tree_data)
        if match:
            sample_tips = sample_tips + match

    # Re-assure pattern search (Sometimes some IDs include others, such as: sample ID LS includes sample ID S...)
    # Check if tip (pattern matched bin base-name), exists in bin dir
    final_tips = list ()
    for tip in sample_tips:
        if (tip+'.fa' in bin_names) or (tip+'_sub.fa' in bin_names):
            final_tips.append(tip)
    final_tips = (',').join(final_tips)


    # Call Rscript to generate sub-trees
    file = os.path.dirname(sys.argv[0])
    curr_dir = os.path.abspath(file)

    subtreeCmd='Rscript '+curr_dir+'/holo-bin_subtre.R --tips '+final_tips+' -in_tree '+tree_path+' -out_tree '+out_tree_path+''
    subprocess.Popen(subtreeCmd,shell=True).wait()







# # Read archea tree
# arc = "cervids_BINF2021/drep_joint/MDR_03-BinPhylogeny/across_CERVIDS/classify/gtdbtk.ar122.classify.tree"
# arc_t = open(arc,'r')
# arc_tree = str()
# for line in arc_t.readlines():
#     arc_tree+=line
# arc_t.close()
#
# # Read bacteria tree
# bac = "cervids_BINF2021/drep_joint/MDR_03-BinPhylogeny/across_CERVIDS/classify/gtdbtk.bac120.classify.tree"
# bac_t = open(bac,'r')
# bac_tree = str()
# for line in bac_t.readlines():
#     bac_tree+=line
# bac_t.close()

# # Find which USER bins are in tree
# bins = list()
# for ID in ID_list:
#     match = re.findall(''+str(ID)+'[0-9]?.{1}[a-z]*_?[0-9]{4}|'+str(ID)+'[0-9]?.{1}[a-z]*_?[0-9]{3}|'+str(ID)+'[0-9]?.{1}[a-z]*_?[0-9]{2}',arc_tree) # creates a match object
#     if match:
#         bins = bins + match
# Check if found bins exist (S and LS bins got confused...)
# final_bins = list()
# real_bins = [os.path.basename(x) for x in glob.glob(dir+'/MDR_01-BinDereplication/across_CERVIDS/dereplicated_genomes/*.fa')]
#
# for bin in bins:
#
#     if bin+'.fa' in real_bins:
#         final_bins.append(bin)
# print(final_bins)
#
