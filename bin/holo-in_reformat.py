#16.04.2020 - Holoflow 0.1.

import subprocess
import argparse
import time
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-r1i', help="read1 input", dest="read1i", required=True)
parser.add_argument('-r2i', help="read2 input", dest="read2i", required=True)
parser.add_argument('-r1o', help="read1 output", dest="read1o", required=True)
parser.add_argument('-r2o', help="read2 output", dest="read2o", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


read1i=args.read1i
read2i=args.read2i
read1o=args.read1o
read2o=args.read2o
ID=args.ID
log=args.log


# Run
if not (os.path.exists(str(read1o))):
    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as log:
        log.write('\t\t'+current_time+'\tInput Files Reformat step - '+ID+'\n')
        log.write('The headers of the .fastq input files are being reformatted.\n\n')


    for i in range(2):
        i+=1
        if i == 1:          # define input output files
            r_i=read1i
            r_o=read1o
        if i == 2:
            r_i=read2i
            r_o=read2o

        with open(str(r_i),'r') as r_input, open(str(r_o), 'w') as r_output:
            n = 1
            read_n=''
            seq1 = ''
            seq2 = ''
            read_id=''
            qual_id=''

            for line in r_input:

                if line.startswith('@'):
                    if seq1:
                        read_n= str(n).zfill(14)
                        read_id = ("@"+str(ID)+"_"+str(read_n)+'.'+str(i)+'\n')
                        r_output.write(read_id+seq1+'\n'+qual_id+seq2+'\n')

                        n += 1
                        seq1=''
                        seq2=''
                        qual_id=''

                    else:
                        pass

                if line.startswith('+'):
                    read_n= str(n).zfill(14)
                    qual_id = ('+\n')

                if seq1 and (not line.startswith('+')):
                    seq2+= line.strip()

                if not (line.startswith('@') or line.startswith('+') or seq2):
                    seq1+= line.strip()


            if seq1:
                read_n= str(n).zfill(14)
                read_id = ("@"+str(ID)+"_"+str(read_n)+'.'+str(i)+'\n')
                r_output.write(read_id+seq1+'\n'+qual_id+seq2+'\n')

                n += 1
                seq1=''
                seq2=''
                qual_id=''

            else:
                pass

if (os.path.isfile(read2o)):
    os.remove(read1i)
    os.remove(read2i)
