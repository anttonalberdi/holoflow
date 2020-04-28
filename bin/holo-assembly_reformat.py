#09.04.2020 - Holoflow 0.1.

import subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-in_a', help="assembly input", dest="in_assembly", required=True)
parser.add_argument('-out_a', help="assembly output", dest="out_assembly", required=True)
parser.add_argument('-st_in', help="stats file input", dest="stats_in", required=True)
parser.add_argument('-st_out', help="stats file output", dest="stats_out", required=True)
args = parser.parse_args()


in_a=args.in_assembly
out_a=args.out_assembly
stats_in=args.stats_in
stats_out=args.stats_out



with open(str(in_a)) as f_input, open(str(out_a), 'w') as f_output:
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


    #Get stats after assembly
    contigs1 = len([1 for line in open(str(in_a)) if line.startswith(">")])

    #Print stats to stats file
    shell('mv '+stats_in+' '+stats_out+'')
    statsfile=open(str(stats_out),"a+")
    statsfile.write("Assembly contigs\t{0} \r\n".format(contigs1))

    #Get stats after assembly reformat
    contigs2 = len([1 for line in open(str(out_a)) if line.startswith(">")])

    #Print stats to stats file
    statsfile.write("Reformated assembly contigs\t{0} \r\n".format(contigs2))
    statsfile.close()
