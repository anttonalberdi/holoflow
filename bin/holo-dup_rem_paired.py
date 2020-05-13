#08.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-o', help="output directory", dest="output_dir", required=True)
parser.add_argument('-sep', help="sep", dest="separator", required=True)
parser.add_argument('-D', help="file to save number and list of dup seqs", dest="file_to_dups")
parser.add_argument('-s', help="by seq", dest="by_seq", required=True)
parser.add_argument('-n', help="by name", dest="by_name", required=True)
parser.add_argument('-i', help="ignore case", dest="ignore", required=True)
args = parser.parse_args()

output_dir=args.output_dir
read1=args.read1
read2=args.read2
separator=args.separator
file_to_dups=args.file_to_dups
by_seq=args.by_seq
by_name=args.by_name
ignore=args.ignore


# Run
if by_seq:
    if (file_to_dups and ignore):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -s -j 28 -o'+ output_dir+' -i -D '+file_to_dups+''

    elif file_to_dups:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -s -j 28 -o'+ output_dir+' -D '+file_to_dups+''

    elif ignore:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -s -j 28 -o'+ output_dir+' -i '

    else:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -s -j 28 -o'+ output_dir+''



if by_name:
    if (file_to_dups and ignore):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -n -j 28 -o'+ output_dir+' -i -D '+file_to_dups+''

    elif file_to_dups:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -n -j 28 -o'+ output_dir+' -D '+file_to_dups+''

    elif ignore:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -n -j 28 -o'+ output_dir+' -i '

    else:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -n -j 28 -o'+ output_dir+''


if not (by_seq or by_name):
    if (file_to_dups and ignore):
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -j 28 -o'+ output_dir+' -i -D '+file_to_dups+''

    elif file_to_dups:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -j 28 -o'+ output_dir+' -D '+file_to_dups+''

    elif ignore:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -j 28 -o'+ output_dir+' -i '

    else:
        seqkitCmd = 'module load tools pigz/2.3.4 seqkit/0.7.1 && paste -d '+separator+' '+read1+' '+read2+' | seqkit rmdup -j 28 -o'+ output_dir+''


subprocess.check_call(seqkitCmd, shell=True)
