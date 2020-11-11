#02.11.2020

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bin_dir', help="drep bin directory", dest="bin_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()



bin_dir=args.bin_dir
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
        logi.write('\t\t'+current_time+'\tBin Annotation step - '+ID+'\n')
        logi.write('\n\n')

    # Get bin names and full paths
    bin_list=glob.glob(str(bin_dir)+"/*.fa")
    for bin in bin_list:
        bin_name=os.path.basename(bin)
        bin=os.path.abspath(bin)

        # Annotation with Prokka
        annCmd='module load tools perl/5.30.2 hmmer/3.2.1 prodigal/2.6.3 tbl2asn/20191211 ncbi-blast/2.8.1+ prokka/1.14.0 && prokka --quiet --cpus '+threads+' --outdir '+out_dir+'/prokka_out --prefix '+bin_name+' '+bin+''
        subprocess.Popen(annCmd, shell=True).wait()

        # Reformat annotations
        functCmd='mkdir '+out_dir+'/bin_funct_annotations && grep product '+out_dir+'/prokka_out/'+bin_name+'/'+bin_name+'.gff > '+out_dir+'/bin_funct_annotations/'+bin_name+'.gff'
        subprocess.check_call(functCmd, shell=True)

        trgenCmd='mkdir '+out_dir+'/bin_translated_genes && cp '+out_dir+'/prokka_out/'+bin_name+'/'+bin_name+'.faa '+out_dir+'/bin_translated_genes'
        subprocess.check_call(trgenCmd, shell=True)

        untrgenCmd='mkdir '+out_dir+'/bin_untranslated_genes && cp '+out_dir+'/prokka_out/'+bin_name+'/'+bin_name+'.ffn '+out_dir+'/bin_untranslated_genes'
        subprocess.check_call(untrgenCmd, shell=True)
