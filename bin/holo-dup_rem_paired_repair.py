#08.04.2020 - Holoflow 0.1

import subprocess
import argparse
import gzip
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-1', help="path1", dest="read1", required=True)
parser.add_argument('-2', help="path2", dest="read2", required=True)
parser.add_argument('-i', help="input_all", dest="input", required=True)
parser.add_argument('-sep', help="sep", dest="separator", required=True)
parser.add_argument('-si', help="stats input file", dest="in_stats", required=True)
parser.add_argument('-so', help="stats output file", dest="out_stats", required=True)
args = parser.parse_args()

input_file=args.input
read1=args.read1
read2=args.read2
separator=args.separator
in_stats=args.in_stats
out_stats=args.out_stats


# Run

# split not dup sequences into reads again: the F-R reads were concatenated in previous rule 
cut1Cmd = 'cut --delimiter='+str(separator)+' -f1 <(zcat '+input_file+') | gzip > '+read1+''
subprocess.Popen(cut1Cmd, shell=True,executable="/bin/bash").wait()
cut2Cmd = 'cut --delimiter='+str(separator)+' -f2  <(zcat '+input_file+') | gzip > '+read2+''
subprocess.Popen(cut2Cmd, shell=True,executable="/bin/bash").wait()
rmCmd = 'rm '+input_file+''
subprocess.check_call(rmCmd, shell=True)


    # Get stats after duplicate removal
mvstatsCmd= 'mv '+in_stats+' '+out_stats+''
subprocess.check_call(mvstatsCmd, shell=True)


reads = 0
bases = 0
with gzip.open(str(read1), 'rt') as read:
  for id in read:
      seq = next(read)
      reads += 1
      bases += len(seq.strip())*2
      next(read)
      next(read)

  #Print stats to stats file
  statsfile=open(str(out_stats),"a+")
  statsfile.write("Dereplicated reads\t{0} ({1} bases)\r\n".format(reads,bases))
  statsfile.close()
