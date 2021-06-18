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
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

a=args.a
bb=args.bb
d=args.d
bt=args.bt
ID=args.ID
log=args.log


# Run
if os.path.exists(bb) and (len(os.listdir(bb)) == 0):
    rmCmd='rm -rf '+bb+''
    subprocess.check_call(rmCmd, shell=True)

if not os.path.exists(bb):

    bin_base = bb+ID+'.vmb'

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tVAMB Binning step - '+ID+'\n')
        logi.write('Individual assembly binning is being done by VAMB. This will sort the contigs into groups,\ncalled bins, which ideally will belong to taxonomically close organisms. This is mainly done\nbased on coverage and tetranucleotide frequencies and differential coverage.\n\n')


# If no bins in directory, then run vamb
    if not glob.glob(str(bb)+"*.fa"):
        vambCmd='module unload gcc && module load tools anaconda3/4.4.0 perl/5.20.2 metabat/2.12.1 && vamb  -o _ --outdir '+bb+' --fasta '+a+' --jgi '+d+' --minfasta 200000'
        subprocess.check_call(vambCmd, shell=True)

            # Modify bin names and create contig to bin table

        binlist=glob.glob(str(bb)+"bins/*.fna")
        n = 0

        for bin in binlist:
            full_bin=os.path.abspath(bin)
            new_bin=bin_base+str(n)+'.fa'
            print(bin)

            renameBinCmd='mv '+full_bin+' '+new_bin+''  # rename to standard
            subprocess.Popen(renameBinCmd, shell=True).wait()
            n +=1

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
            with open(bin_base+'_checked_bins','w+') as check:
                check.write('True Vamb vmb')

        else:
            with open(bin_base+'_checked_bins','w+') as check:
                check.write('False Vamb vmb')


        os.rmdir(bb+'bins')
