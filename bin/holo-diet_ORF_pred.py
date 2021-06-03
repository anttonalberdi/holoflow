#06.05.2021 - Holoflow 0.1.
import subprocess
import argparse
import os
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-out_dir', help="out_dir", dest="out_dir", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


a=args.a
out_dir=args.out_dir
t=args.threads
ID=args.ID
log=args.log



# Run
# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\tHOLOFLOW\tMETAGENOMICS\n\t\t'+current_time+'\t - '+ID+'\n')
    logi.write('Genes are being predicted by Prodigal 2.6.3.\n\n')


# Generate .faa and .fna outputs
out_coords = out_dir+'/'+ID+'.coords.gff'
ptranslations = out_dir+'/'+ID+'.ptranslations.faa'
nsequences = out_dir+'/'+ID+'.predORFs.fna'

if not os.path.isfile(ptranslations):
    prodigalCmd='module unload gcc && module load tools prodigal/2.6.3 && prodigal -i '+a+' -o '+out_coords+' -a '+ptranslations+' -p meta -f gff -d '+nsequences+''
    subprocess.check_call(prodigalCmd, shell=True)
