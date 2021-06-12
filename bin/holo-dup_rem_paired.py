#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import os


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-o ', help="output directory", dest="output", required=True)
parser.add_argument('-sep', help="sep", dest="separator", required=True)
parser.add_argument('-D', help="file to save number and list of dup seqs", dest="file_to_dups",required=True)
parser.add_argument('-s', help="by seq", dest="by_seq", required=True)
parser.add_argument('-n', help="by name", dest="by_name", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-i', help="ignore case", dest="ignore", required=True)
args = parser.parse_args()

output=args.output
read1=args.read1
read2=args.read2
separator=args.separator
file_to_dups=args.file_to_dups
by_seq=args.by_seq
by_name=args.by_name
ID=args.ID
log=args.log
ignore=args.ignore


# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as log:
    log.write('\t\t'+current_time+'\tDuplicates Removal step - '+ID+'\n')
    log.write('Duplicate sequences are being removed.\n\n')

# compressed input and outputs
# all different conditions for different variables in config that can be used, modified or not used at all. Not very optimal
if by_seq == 'True':
    if (not file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -s -i -D '+file_to_dups+' -o '+ output+''

    elif (not file_to_dups == 'False') and (ignore == 'False'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -s -D '+file_to_dups+' -o '+ output+''

    elif (file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -s -i -o '+ output+''

    else:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -s -o '+ output+''



if by_name == 'True':
    if (not file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -n -i -D '+file_to_dups+' -o '+output+''

    elif (not file_to_dups == 'False') and (ignore == 'False'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -n -D '+file_to_dups+' -o '+output+''

    elif (file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -n -i -o '+output+''

    else:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' <(zcat '+read1+') <(zcat '+read2+') | seqkit -j 40 rmdup -n -o '+output+''

subprocess.check_call(seqkitCmd, shell=True,executable="/bin/bash")
