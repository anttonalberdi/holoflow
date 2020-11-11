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
            bin_name=bin
            bin=os.path.abspath(bin)

        # Annotation with Prokka
        ######### DEPENDENCIES module load perl/5.30.2 hmmer/3.2.1 TEST MORE
        annCmd='prokka --quiet --cpus '+threads+' --outdir '+out_dir+' --prefix '+bin_name+' '+bin+''
        subprocess.check_call(annCmd, shell=True)






for i in $(ls ${bins}); do
        bin_name=${i%.*}
        bin_file=${bins}/$i
	echo "${SOFT}/shorten_contig_names.py $bin_file > ${out}/tmp_bin.fa"
	${SOFT}/shorten_contig_names.py $bin_file > ${out}/tmp_bin.fa
	if [[ $? -ne 0 ]]; then error "Could not process/shorten the contig names of ${bin_file}. Exiting..."; fi
        comm "NOW ANNOTATING ${bin_name}"

        cmd="prokka --quiet --cpus $threads --outdir ${out}/prokka_out/$bin_name --prefix $bin_name ${out}/tmp_bin.fa"
	echo $cmd
	$cmd

	if [[ $? -ne 0 ]]; then warning "Something possibly went wrong with annotating ${bin_name}. Proceeding anyways"; fi
        if [[ ! -s ${out}/prokka_out/${bin_name}/${bin_name}.gff ]]; then error "Something went wrong with annotating ${bin_name}. Exiting..."; fi
	rm ${out}/tmp_bin.fa
done



if [[ $(ls ${out}/prokka_out/) -lt 1 ]]; then error "Something went wrong with running prokka on all the bins! Exiting..."; fi

comm "PROKKA finished annotating all the bins!"





    if (os.path.exists(str(''+out_dir+'/final_bins_Info.csv'))):
        drepbinsCmd=''
        subprocess.check_call(drepbinsCmd, shell=True)
