#20.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import time
import re

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-d', help="depth file", dest="d", required=True)
parser.add_argument('-bb', help="bin base ID", dest="bb", required=True)
parser.add_argument('-bt', help="bin table output", dest="bt", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

a=args.a
d=args.d
bb=args.bb
bt=args.bt
t=args.t
ID=args.ID
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\t\t'+current_time+'\tMaxbin Binning step - '+ID+'\n')
    logi.write('Individual assembly binning is being done by MAXBIN. This will sort the contigs into groups,\ncalled bins, which ideally will belong to taxonomically close organisms. This is mainly done\nbased on coverage and tetranucleotide frequencies.\n\n')




if not glob.glob(str(bb)+"*.fa"):
    maxbinCmd='module unload gcc && module load tools perl/5.20.2 maxbin/2.2.7 fraggenescan/1.31 && run_MaxBin.pl -contig '+a+' -abund '+d+' -out '+bb+' -thread '+t+''
    subprocess.check_call(maxbinCmd, shell=True)

        # Modify bin names and create contig to bin table
    renamebinsCmd='binlist=$(ls '+bb+'*.fasta | sed "s/.*mxb\.//" | sed "s/\.fasta//") && for bin in $binlist; do bin2=$((10#$bin)) ; mv '+bb+'.${bin}.fasta '+bb+'${bin2}.fa; done'
    subprocess.Popen(renamebinsCmd, shell=True).wait()


        #Fill contig to bin table
    binlist=glob.glob(str(bb)+"*.fa")
    bintable = open(str(bt),"a+")

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
            check.write('True maxbin mxb')

    else:
        with open(bb+'_checked_bins','w+') as check:
            check.write('False maxbin mxb')
