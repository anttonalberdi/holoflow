#20.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-coa', help="coassembly TRUE or FALSE", dest="coa", required=True)
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-d', help="depth file", dest="d", required=True)
parser.add_argument('-bb', help="bin base ID", dest="bb", required=True)
parser.add_argument('-bt', help="bin table output", dest="bt", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-l', help="minimum contig length", dest="l", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

coa=args.coa
a=args.a
d=args.d
bb=args.bb
bt=args.bt
t=args.t
l=args.l
log=args.log

# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tConcoct Binning step\n')
    log.write('Coassembly binning is being done by CONCOCT. (((MERGE SAMPLES))) This will sort the contigs into groups,\ncalled bins, which ideally will belong to taxonomically close organisms. This is mainly done\nbased on coverage and tetranucleotide frequencies.\n\n')


if coa: # default set to FALSE in configfile    # first bin 0 --> to +1 
    if not glob.glob(str(bb)+"*.fa"):
        concoctCmd='concoct --coverage_file '+d+' --composition_file '+a+' -b '+bb+' -l '+int(l)+''
        subprocess.check_call(concoctCmd, shell=True)

            #Create contig to bin table
        bintable = open(str(bt),"a+")
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
