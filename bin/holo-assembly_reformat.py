#09.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="assembly", required=True)
parser.add_argument('-o', help="output directory", dest="output", required=True)
args = parser.parse_args()


out=args.out
assembly=args.assembly


# Reformat contig names and filter by contig length
with open(str(assembly)) as f_input, open(str(output), 'w') as f_output:
    seq = ''
    contig_n = 0

    for line in f_input:
        if line.startswith('>'):

            if seq:
                if len(seq) > 1000:
                    contig_n += 1
                    contig_id = (">C_"+str(contig_n))
                    seq += ('\n')

                    f_output.write(contig_id + '\n' + seq)
                    seq = ''

                else:
                    seq = ''
        else:
            seq += line.strip()

    if seq:
        if len(seq) > 1000:
            contig_n += 1
            contig_id = (">C_"+str(contig_n))
            seq += ('\n')
            f_output.write(contig_id + '\n' + seq)

        else:
            pass
