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
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMetabat Binning step - '+ID+'\n')
    log.write('Individual assembly binning is being done by METABAT. This will sort the contigs into groups,\ncalled bins, which ideally will belong to taxonomically close organisms. This is mainly done\nbased on coverage and tetranucleotide frequencies.\n\n')



if not glob.glob(str(bb)+"*.fa"):
    metabatCmd='module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && metabat2 -i '+a+' -a '+d+' -o '+bb+' -m 1500 -t '+t+''
    subprocess.Popen(metabatCmd, shell=True).wait()

        #Fill contig to bin table
    binlist=glob.glob(str(bb)+"*.fa")
    bintable = open(str(bt),"a+")

    for bin in binlist:
        full_bin=os.path.abspath(bin)
        new_bin=full_bin.replace("mtb.","mtb")

        renameBinCmd='mv '+full_bin+' '+new_bin+''
        subprocess.check_call(renameBinCmd, shell=True)

    binlist=glob.glob(str(bb)+"*.fa")
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
            check.write('True metabat mtb')

    else:
        with open(bb+'_checked_bins','w+') as check:
            check.write('False metabat mtb')
