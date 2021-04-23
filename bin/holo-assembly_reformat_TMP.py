#09.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import os


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-in_a', help="assembly input", dest="in_assembly", required=True)
parser.add_argument('-out_a', help="assembly output", dest="out_assembly", required=True)
parser.add_argument('-st_in', help="stats file input", dest="stats_in", required=True)
parser.add_argument('-st_out', help="out directory", dest="out", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-min_cl', help="minimum contig length", dest="min_cl", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


in_a=args.in_assembly
out_a=args.out_assembly
stats_in=args.stats_in
ID=args.ID
min_cl=args.min_cl
out=args.out
log=args.log


# Run
if os.path.exists(str(out_a)):
    pass

if not os.path.exists(str(out_a)):

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as log:
        log.write('\t\t'+current_time+'\tAssembly Reformat step - '+ID+'\n')
        log.write('The generated assembly file in the previous step is being reformatted: Those contigs less than '+min_cl+'\nbase pairs long are being removed and the IDs of the remaining ones are being modified.\n\n')

    #unzip temp assembly
    unzCmd='gunzip '+in_a+''
    subprocess.Popen(unzCmd,shell=True).wait()
    in_a = in_a.replace('.gz','')
    out_a = out_a.replace('.gz','')


    with open(str(in_a)) as f_input, open(str(out_a), 'w') as f_output:
        seq = ''
        # create list with six-digit numbers: 000001 -> 100000
        # to re-enumerate the contigs
        contig_n = (["%06d" % x for x in range(1000000)])
        n = 0

        # the assembly has two lines per contig : > ID and sequence
        for line in f_input:
            if line.startswith('>'):
                # If the line corresponds to the ID, create new ID with 6-digit numeration + group ID
                # for the PREVIOUS contig. This loop only stores in variables the SEQUENCES, so for
                # every sequence, a new contig ID is generated

                if seq:
                    # Carry on only if the sequence paired with this ID is longer than the minimum contig length
                    # provided by the user - default 1500bp, otherwise continue and omit this contig
                    if len(seq) > int(min_cl):
                        n += 1
                        contig_id = (">"+str(ID)+"_"+str(contig_n[n]))
                        # add new line after sequence
                        seq += ('\n')
                        # Write to new assembly reformatted file
                        f_output.write(contig_id + '\n' + seq)
                        # un-define sequence, and continue to next
                        seq = ''

                    else:
                        seq = ''
            else:
                seq += line.strip()
        # Last line - the loop has finished but the last contig has not yet been reformatted + written
        if seq:
            if len(seq) > int(min_cl):
                n += 1
                contig_id = (">"+str(ID)+"_"+str(contig_n[n]))
                seq += ('\n')
                f_output.write(contig_id + '\n' + seq)

            else:
                pass


    #Get stats after assembly
    contigs1 = len([1 for line in open(str(in_a)) if line.startswith(">")])

    #Print stats to stats file

    statsfile=open(str(stats_in),"a+")
    statsfile.write("Assembly contigs\t"+str(contigs1)+" \r\n")

    #Get stats after assembly reformat
    contigs2 = len([1 for line in open(str(out_a)) if line.startswith(">")])

    #Print stats to stats file
    statsfile.write("Reformated assembly contigs\t"+str(contigs2)+" \r\n")
    statsfile.close()

    statsCmd='mv '+stats_in+' '+out+''
    subprocess.check_call(statsCmd, shell=True)

    # gzip outputs
    gzCmd='gzip '+in_a+' '+out_a+''
    subprocess.Popen(gzCmd,shell=True).wait()


else:
    pass
