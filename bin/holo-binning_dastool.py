#27.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-bt_mtb', help="metabat bin table", dest="bt_mtb", required=True)
parser.add_argument('-bt_mxb', help="maxbin bin table", dest="bt_mxb", required=True)
parser.add_argument('-fbt', help="temp joined bin table", dest="fbt", required=True)
parser.add_argument('-p', help="prodigal predicted proteins", dest="p", required=True)
parser.add_argument('-o', help="output main dir", dest="o", required=True)
parser.add_argument('-bin_o', help="bin final dir", dest="bin_o", required=True)
parser.add_argument('-se', help="search engine", dest="se", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-db', help="dastool database directory", dest="db", required=True)
args = parser.parse_args()


a=args.a
bt_mtb=args.bt_mtb
bt_mxb=args.bt_mxb
fbt=args.fbt
p=args.p
o=args.o
bin_o=args.bin_o
se=args.se
t=args.t
db=args.db



# Run

dastoolDependencies='module load tools gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ prodigal/2.6.3 das_tool/1.1.1 diamond/0.9.24 usearch/11.0.667'
# with open(str(bt_mxb), 'r') as mxb, open(str(bt_mtb), 'r') as mtb:
#     bincontig_tables = mtb.readlines()
#     bincontig_tables = bincontig_tables.append(mxb.readlines())
#
#
# # Reading bin tables to temp bin table
# # with open(str(bt_mxb), 'r') as mxb:
# #     r1 = mxb.read()
# #
# # with open(str(bt_mtb), 'r') as mtb:
# #     r2 = mtb.read()
# #
# # r1 += "\n"
# # r1 += r2
# #
# # createtempCmd='touch '+tbt+''
# # subprocess.check_call(createtempCmd, shell=True)
# #
# # with open (str(tbt), 'w+') as tbt_w:
# #     tbt_w.write(r1)

bincontig_tables = ",".join(glob.glob(os.path.join(fbt, '_*.txt')))

dastoolCmd=''+dastoolDependencies+' DAS_Tool -i '+bincontig_tables+' -c '+a+' -o '+o+' --proteins '+p+' -l maxbin,metabat --search_engine '+se+' -t '+t+' --db_directory '+db+' --write_bins 1'
subprocess.check_call(dastoolCmd, shell=True)
#removetempCmd='rm '+tbt+''
#subprocess.check_call(removetempCmd, shell=True)


# Move definitive bins to final directory
binfiles = glob.glob(os.path.join(str(o),'*.fa'))
for b in binfiles:
    shutil.move(b, str(bin_o))
