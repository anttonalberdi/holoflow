#22.11.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import glob
import numpy as np
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bam_dir', help="input bam from mapped MAGs to .fastq directory", dest="bam_dir", required=True)
parser.add_argument('-mag_dir', help="originally dereplicated mags", dest="mag_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

bam_dir=args.bam_dir
mag_dir=args.mag_dir
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
        logi.write('\t\t'+current_time+'\tMAG Coverage step - '+ID+'\n')
        logi.write('\n\n')

    # # Extract MAGs coverage from bam files - BY CONTIG
    #     # CONTIGS X SAMPLES
    depth_contig=out_dir+'/'+ID+'.coverage_byContig.txt'
    getcoverageCmd='module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+depth_contig+' '+str(bam_dir)+'/*.bam'
    subprocess.check_call(getcoverageCmd, shell=True)

    # Generate aggregated coverage table  - BY MAG
        # MAGS X SAMPLES
    depth_mag=out_dir+'/'+ID+'.coverage_byMAG.txt'
    coverage_data=list()

    with open(depth_mag, 'w+') as cov_mag:

        # Start MAG table with same line as depth_mag
        cov_contig = open(depth_contig,'r')
        first_dcontig = cov_contig.readline()
        first_dcontig = first_dcontig.replace('contig','MAG')
        cov_mag.write(first_dcontig.strip()+'\n')
        cov_contig.close()

        # Prepare mag data and ID
        mag_list=glob.glob(str(mag_dir)+'/*.fa')
        for mag in mag_list:
            mag_id=''
            cov_data_tomag=''
            mag_id=os.path.basename(mag)
            mag_id=mag_id.replace('.fa','')
            if '.contigs' in mag_id:
                mag_id=mag_id.replace('.contigs','')

            tmp_MAGcoverage=out_dir+'/'+ID+'.'+mag_id+'_MAGcoverage.txt'

            grepCmd='grep '+mag_id+' '+depth_contig+' > '+tmp_MAGcoverage+''
            subprocess.Popen(grepCmd, shell=True).wait()

            # Sum coverage and length stats for contigs in same mag, write
            cov_data_id=np.genfromtxt(tmp_MAGcoverage,delimiter='\t')
            cov_data_id=np.array(cov_data_id)
            cov_data = np.delete(cov_data_id, 0, 1) # remove contig ID column

            # Sum coverage and length for all contigs in mag
            cov_data=cov_data.astype(np.float)
            cov_data=np.sum(cov_data,axis=0)
            cov_data=cov_data.round(decimals=4)
            cov_data=cov_data.tolist()

            # Write coverage for given MAG
            for num in cov_data:
                cov_data_tomag+=str(num)+'\t'

            cov_mag.write(mag_id+'\t'+str(cov_data_tomag)+'\n')
            os.remove(tmp_MAGcoverage)
