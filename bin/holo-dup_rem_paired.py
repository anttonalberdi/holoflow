#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-o ', help="output directory", dest="output_dir", required=True)
parser.add_argument('-sep', help="sep", dest="separator", required=True)
parser.add_argument('-D', help="file to save number and list of dup seqs", dest="file_to_dups")
parser.add_argument('-s', help="by seq", dest="by_seq", required=True)
parser.add_argument('-n', help="by name", dest="by_name", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-i', help="ignore case", dest="ignore", required=True)
#parser.add_argument('--foo', action='store_true') - would be the optimal option if not Snakemake

args = parser.parse_args()

output_dir=args.output_dir
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
    log.write('\t\t'+current_time+'\tDuplicates Removal step - ID '+ID+'\n')
    log.write('Duplicate sequences are being removed.\n\n')



if by_seq == 'True':

    if (not file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -s -i -D '+file_to_dups+' -o '+ output_dir+''

    elif (not file_to_dups == 'False') and (ignore == 'False'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -s -D '+file_to_dups+' -o '+ output_dir+''

    elif (file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -s -i -o '+ output_dir+''

    else:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -s -o '+ output_dir+''



if by_name == 'True':
    if (not file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -n -i -D '+file_to_dups+' -o '+ output_dir+''

    elif (not file_to_dups == 'False') and (ignore == 'False'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -n -D '+file_to_dups+' -o '+ output_dir+''

    elif (file_to_dups == 'False') and (ignore == 'True'):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -n -i -o '+ output_dir+''

    else:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -n -o '+ output_dir+''

print(seqkitCmd)
subprocess.check_call(seqkitCmd, shell=True)


# if not (by_seq or by_name):
#     if (file_to_dups and ignore):
#         seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -o '+ output_dir+' -i -D '+file_to_dups+''
#
#     if (not ignore) and file_to_dups:
#         seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -o '+ output_dir+' -D '+file_to_dups+''
#
#     if (not file_to_dups) and ignore:
#         seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -o '+ output_dir+' -i '
#
#     else:
#         seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit -j 40 rmdup -o '+ output_dir+''
#
