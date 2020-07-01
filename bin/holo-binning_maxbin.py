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
parser.add_argument('-sample', help="sample", dest="sample", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

a=args.a
d=args.d
bb=args.bb
bt=args.bt
t=args.t
sample=args.sample
log=args.log


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tMaxbin Binning step - Sample '+sample+'\n')
    log.write('Individual assembly binning is being done by MAXBIN. This will sort the contigs into groups,\ncalled bins, which ideally will belong to taxonomically close organisms. This is mainly done\nbased on coverage and tetranucleotide frequencies.\n\n')




if not glob.glob(str(bb)+"*.fasta"):
    try:
        maxbinCmd='module unload gcc && module load tools perl/5.20.2 maxbin/2.2.7 fraggenescan/1.31 && run_MaxBin.pl -contig '+a+' -abund '+d+' -out '+bb+' -thread '+t+''
        subprocess.check_call(maxbinCmd, shell=True)

            #Create contig to bin table
        bintable = open(str(bt),"a+")
        binlist=glob.glob(str(bb)+"*.fasta")


        for bin in binlist:
            binname = os.path.splitext(os.path.basename(bin))[0]+''
            with open(bin, 'r') as binfile:
               for line in binfile:
                    if line.startswith('>'):
                        contig = line.strip()
                        contig = contig.replace(">", "")
                        bintable.write("{0}\t{1}\r\n".format(contig,binname))
        bintable.close()

    except:
        # Write to log
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        with open(str(log),'a+') as log:
            log.write(''+current_time+' - Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop.\n\n')
        pass
