#20.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-d', help="depth file", dest="d", required=True)
parser.add_argument('-bb', help="bin base ID", dest="bb", required=True)
parser.add_argument('-bt', help="bin table output", dest="bt", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
args = parser.parse_args()

a=args.a
d=args.d
bb=args.bb
bt=args.bt
t=args.t



if not glob.glob(str(bb)+"*.fasta"):
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
