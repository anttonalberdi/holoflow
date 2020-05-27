#27.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-dep', help="dastool dependencies", dest="dep", required=True)
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-bt_mtb', help="metabat bin table", dest="bt_mtb", required=True)
parser.add_argument('-bt_mxb', help="maxbin bin table", dest="bt_mxb", required=True)
parser.add_argument('-p', help="prodigal predicted proteins", dest="p", required=True)
parser.add_argument('-o', help="output main dir", dest="o", required=True)
parser.add_argument('-bin_o', help="bin final dir", dest="bin_o", required=True)
parser.add_argument('-se', help="search engine", dest="se", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-db', help="dastool database directory", dest="db", required=True)
args = parser.parse_args()


dep=args.dep
a=args.a
bt_mtb=args.bt_mtb
bt_mxb=args.bt_mxb
p=args.p
o=args.o
bin_o=args.bin_o
se=args.se
t=args.t
db=args.db



# Run

bincontig_tables=",".join(glob.glob(str(bt_mxb),str(bt_mtb)))
dastoolCmd=''+de+' && DAS_Tool -i '+bincontig_tables+' -c '+a+' -o '+o+' --proteins '+p+' -l maxbin,metabat --search_engine '+se+' -t '+t+' --db_directory '+db+' --write_bins 1'
subprocess.check_call(dastoolCmd, shell=True)


# Move definitive bins to final directory

binfiles = glob.glob(os.path.join(str(o),'*.fa'))
for b in binfiles:
    shutil.move(b, str(bin_o))
