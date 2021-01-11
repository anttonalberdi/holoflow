#27.05.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import sys
import glob
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-a', help="assembly file", dest="a", required=True)
parser.add_argument('-cb', help="checked bins", dest="check_b")
parser.add_argument('-bt_mtb', help="metabat bin table", dest="bt_mtb", required=True)
parser.add_argument('-bt_mxb', help="maxbin bin table", dest="bt_mxb", required=True)
parser.add_argument('--bt_cct', help="concoct bin table", dest="bt_cct")
#parser.add_argument('-p', help="prodigal predicted proteins", dest="p", required=True)
parser.add_argument('-o', help="output main dir", dest="o", required=True)
parser.add_argument('-se', help="search engine", dest="se", required=True)
parser.add_argument('-t', help="threads", dest="t", required=True)
parser.add_argument('-db', help="dastool database directory", dest="db", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()

a=args.a
bt_mtb=args.bt_mtb
bt_mxb=args.bt_mxb
#p=args.p
o=args.o
se=args.se
t=args.t
db=args.db
ID=args.ID
log=args.log



# Run

# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\t\t'+current_time+'\tDASTool Bin Refinement step - '+ID+'\n')
    logi.write('The binning results from MaxBin and Metabat2 are integrated by DASTool to produce one only non-redundant\nset of bins between them.\n\n')

if args.check_b: # means all binners have bins, either duplicated or own
    bin_dir=os.path.dirname(bt_mtb)
    rmCmd='rm -rf '+args.check_b+' '+bin_dir+'/*remove'
    subprocess.check_call(rmCmd,shell=True)

    # Coassembly
    if args.bt_cct:
        bt_cct=args.bt_cct

        dastoolDependencies='module unload maxbin/2.2.7 fraggenescan/1.31 perl/5.20.2 && module load tools gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ prodigal/2.6.3 das_tool/1.1.1 diamond/0.9.24 usearch/11.0.667'
        dastoolCmd=''+dastoolDependencies+' && DAS_Tool -i '+bt_cct+','+bt_mxb+','+bt_mtb+' -c '+a+' -o '+o+' -l concoct,maxbin,metabat --search_engine '+se+' -t '+t+' --db_directory '+db+' --write_bins 1'
        #dastoolCmd=''+dastoolDependencies+' && DAS_Tool -i '+bt_cct+','+bt_mxb+','+bt_mtb+' -c '+a+' -o '+o+' --proteins '+p+' -l concoct,maxbin,metabat --search_engine '+se+' -t '+t+' --db_directory '+db+' --write_bins 1'
        subprocess.check_call(dastoolCmd, shell=True)


        # Move definitive bins to final directory
        # Remove '.contigs' from bin ID, which was added by DASTool
        ori_dir=o+"_DASTool_bins"
        bins=glob.glob(ori_dir+"/*.fa")

        for bin in bins:
          new_bin=bin.replace('.contigs','')

          if not (new_bin == bin):
            renameCmd='mv '+bin+' '+new_bin+''
            subprocess.check_call(renameCmd,shell=True)

        # Move definitive bins to final directory
        bins=glob.glob(o+"_DASTool_bins/*.fa")

        for bin in bins:
            mvCmd='cd '+ori_dir+'/.. && mv '+bin+' . && rm -rf '+ori_dir+''
            subprocess.check_call(mvCmd,shell=True)


        if os.path.exists(str(o+'_maxbin.eval')):
            # Add relevant info to log
            with open(str(log),'a+') as logf:

                logf.write('\t\tDASTool MaxBin bins evaluation - ID '+ID+'\n\n')
                with open(str(o+'_maxbin.eval'),'r') as mxb_eval:
                    logf.write(''+mxb_eval.read()+'\n\n\n')

                logf.write('\t\tDASTool Metabat2 bins evaluation - ID '+ID+'\n\n')
                with open(str(o+'_metabat.eval'),'r') as mtb_eval:
                    logf.write(''+mtb_eval.read()+'\n\n\n')

                logf.write('\t\tDASTool Concoct bins evaluation - ID '+ID+'\n\n')
                with open(str(o+'_concoct.eval'),'r') as cct_eval:
                    logf.write(''+cct_eval.read()+'\n\n\n')

                if os.path.exists(str(o+'_DASTool_summary.txt')):
                    logf.write('\t\tDASTool Bin Merging Summary - ID '+ID+'\n\n')
                    with open(str(o+'_DASTool_summary.txt'),'r') as summary:
                        logf.write(''+summary.read()+'\n\n\n\n')
                else:
                    pass


    else: # Individual assembly and binning - only maxbin and metabat

        dastoolDependencies='module unload maxbin/2.2.7 fraggenescan/1.31 perl/5.20.2 && module load tools gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ prodigal/2.6.3 das_tool/1.1.1 diamond/0.9.24 usearch/11.0.667'
        dastoolCmd=''+dastoolDependencies+' && DAS_Tool -i '+bt_mxb+','+bt_mtb+' -c '+a+' -o '+o+' --proteins '+p+' -l maxbin,metabat --search_engine '+se+' -t '+t+' --db_directory '+db+' --write_bins 1'
        subprocess.check_call(dastoolCmd, shell=True)


        # Remove '.contigs' from bin ID, which was added by DASTool
        ori_dir=o+"_DASTool_bins"
        bins=glob.glob(ori_dir+"/*.fa")

        for bin in bins:
          new_bin=bin.replace('.contigs','')

          if not (new_bin == bin):
            renameCmd='mv '+bin+' '+new_bin+''
            subprocess.check_call(renameCmd,shell=True)

        # Move definitive bins to final directory and rest to sub-dir
        bins=glob.glob(o+"_DASTool_bins/*.fa")

        for bin in bins:
            # bins in DASTool bins and rest of files in DASTool files && bins out to main dir, remove DASTool bins dir
            mvCmd='mv '+o+'_DASTool_summary.txt '+o+'_DASTool_bins && mkdir '+ori_dir+'/DASTool_files && find '+ori_dir+' -maxdepth 1 -type f | xargs -I {} cp {} '+ori_dir+'/DASTool_files && mv '+o+'_DASTool_bins/* '+ori_dir+' && rm -rf '+o+'_DASTool_bins'
            subprocess.check_call(mvCmd,shell=True)


        # Write to log
        if os.path.exists(str(o+'_maxbin.eval')):
            # Add relevant info to log
            with open(str(log),'a+') as logf:

                logf.write('\t\tDASTool MaxBin bins evaluation - ID '+ID+'\n\n')
                with open(str(o+'_maxbin.eval'),'r') as mxb_eval:
                    logf.write(''+mxb_eval.read()+'\n\n\n')

                logf.write('\t\tDASTool Metabat2 bins evaluation - ID '+ID+'\n\n')
                with open(str(o+'_metabat.eval'),'r') as mtb_eval:
                    logf.write(''+mtb_eval.read()+'\n\n\n')

                if os.path.exists(str(o+'_DASTool_summary.txt')):
                    logf.write('\t\tDASTool Bin Merging Summary - ID '+ID+'\n\n')
                    with open(str(o+'_DASTool_summary.txt'),'r') as summary:
                        logf.write(''+summary.read()+'\n\n\n\n')
                else:
                    pass

else: # No binners had bins
    sys.exit()
