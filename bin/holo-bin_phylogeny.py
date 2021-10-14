#11.11.2020

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-genome_dir', help="genomes directory", dest="gen_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


gen_dir=args.gen_dir
gen_dir=str(gen_dir+"/dereplicated_genomes")
out_dir=args.out_dir
ID=args.ID
log=args.log
threads=args.threads

# Run
if not (os.path.exists(str(out_dir))):
    os.mkdir(str(out_dir))

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tTaxonomic Classification step - '+ID+'\n')
        logi.write('GTDB-Tk is assigning objective taxonomic classifications to baterial genomes based on the Genome Database Taxonomy GTDB.\nThe taxonomic classifications can be found in the .summary.tsv file.\n\n')

    # Call gtdbtk
    gtdbtkCmd='module load tools anaconda3/4.4.0 anaconda2/4.4.0 gtdbtk/1.6.0 prodigal/2.6.3 hmmer/3.2.1 pplacer/1.1.alpha19 fastani/1.1 && \
    gtdbtk classify_wf --genome_dir '+gen_dir+' --extension "fa" --out_dir '+out_dir+' --cpus '+threads+' --scratch_dir '+out_dir+''
    subprocess.Popen(gtdbtkCmd,shell=True).wait()
