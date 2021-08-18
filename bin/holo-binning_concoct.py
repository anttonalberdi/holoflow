#!/usr/bin/env python3

#20.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-d', help="depth file", dest="d", required=True)
parser.add_argument('-bb', help="bin base ID", dest="bb", required=True)
parser.add_argument('-bt', help="bin table output", dest="bt", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-l', help="minimum contig length", dest="l", required=True)
parser.add_argument('-r', help="minimum contig length", dest="r", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
d=args.d
bb=args.bb
bt=args.bt
t=args.t
l=args.l
r=args.r
ID=args.ID
log=args.log

# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tConcoct Binning step\n')
    log.write('Coassembly binning is being done by CONCOCT. This will sort the contigs into groups,\ncalled bins, which ideally will belong to taxonomically close organisms. This is mainly done\nbased on coverage and tetranucleotide frequencies.\n\n')

output_path=bb.replace('/'+ID+'.cct','')

# If the output directory does not contain bins, continue
if not glob.glob(output_path+"/*.fa"):
    # Three-step concoct run
    if not os.path.isfile(''+bb+'_PCA_components_data_gt1500.csv'):
        concoct1Cmd='module load tools && concoct --coverage_file '+d+' --no_original_data --composition_file '+a+' -b '+bb+' -l '+l+' -t '+t+' -r '+r+' '
        subprocess.Popen(concoct1Cmd, shell=True).wait()
    else:
        pass

    if not os.path.isfile(''+bb+'_clustering_merged.csv'):
        concoct2Cmd='merge_cutup_clustering.py '+bb+'_clustering_gt1500.csv > '+bb+'_clustering_merged.csv  && mv '+bb+'_clustering_merged.csv? '+bb+'_clustering_merged.csv' # The script creates ? in the end of the name file: Sounds like you script uses \r\n as line endings, this is typical DOS style line endings. Unix like systems uses \n.
        subprocess.Popen(concoct2Cmd, shell=True).wait()
    else:
        pass

    concoct3Cmd='extract_fasta_bins.py '+a+' '+bb+'_clustering_merged.csv --output_path '+output_path+''
    subprocess.Popen(concoct3Cmd, shell=True).wait()


        #Create contig to bin table
    bintable = open(str(bt),"a+")

    # Rename bins
    binlist=glob.glob(output_path+"/*.fa")

    # Rename bins to standard 
    for bin in binlist:
        full_bin=os.path.abspath(bin)
        base_bin=os.path.basename(bin)
        new_bin=bb+base_bin

        renameBinCmd='mv '+full_bin+' '+new_bin+''
        subprocess.check_call(renameBinCmd, shell=True)


    binlist=glob.glob(bb+'*.fa')

    for bin in binlist:
        binname = os.path.splitext(os.path.basename(bin))[0]+''
        with open(bin, 'r') as binfile:
           for line in binfile:
                if line.startswith('>'):
                    contig = line.strip()
                    contig = contig.replace(">", "")
                    bintable.write("{0}\t{1}\r\n".format(contig,binname))
    bintable.close()


# check
    if binlist: # if bin list not empty, which means bin table exists
        with open(bb+'_checked_bins','w+') as check:
            check.write('True concoct cct')

    else:
        with open(bb+'_checked_bins','w+') as check:
            check.write('False concoct cct')
