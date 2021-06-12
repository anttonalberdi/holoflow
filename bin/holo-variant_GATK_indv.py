## 11.01.20 - Holoflow 0.1

import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bam_dir', help="bam files directory", dest="bam_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-ref_g', help="reference genome", dest="ref_g", required=True)
parser.add_argument('-min_pruning', help="minimum pruning", dest="min_pruning", required=True)
parser.add_argument('-min_dangling', help="minimum dangling", dest="min_dangling", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()


bam_dir=args.bam_dir
out_dir=args.out_dir
ref_g=args.ref_g
min_pruning=args.min_pruning
min_dangling=args.min_dangling
chr_list=args.chr_list
ID=args.ID
log=args.log
threads=args.threads

## Run
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tVariant calling with GATK step - '+ID+'\n')
        logi.write(' \n\n')

    # Get chromosomes list
    chromosome_list = list()
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            chromosome_list.append(chr.strip())

    # Generate empty sample map files per each chromosome
    for CHR in chromosome_list:
        sample_map_file = out_dir+'/sample_map.'+CHR+'.txt'
        os.mknod(sample_map_file)

    # Create .dict for reference genome
    ref_g_base = ref_g.replace('.fna','')
    if not os.path.isfile(ref_g_base+'.dict'):
        dictCmd = 'PICARD="/services/tools/picard-tools/2.20.2/picard.jar" && java -jar $PICARD CreateSequenceDictionary R='+ref_g+' O='+ref_g_base+'.dict'
        subprocess.Popen(dictCmd,shell=True).wait()


    # Generate bam files' paths list & index
    bam_list = glob.glob(bam_dir+'/*.bam')

    for bam in bam_list:
        bam_ID = bam.replace(bam_dir+'/','')
        bam_ID = bam_ID.replace('.bam','')
        if '_ref' in bam_ID:
            bam_ID = bam_ID.replace('_ref','')

        # Index bam with samtools
        if not os.path.isfile(bam+'.bai'):
            idxCmd = 'module load tools samtools/1.11 && samtools index '+bam+''
            subprocess.Popen(idxCmd,shell=True).wait()


        for CHR in chromosome_list:
            out_haplo = out_dir+'/'+bam_ID+'_'+CHR+'.raw.g.vcf.gz'

            if not (min_dangling == 'False'):

                if not (min_pruning == 'False'):
                    haploCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk HaplotypeCaller --java-options "-Xmx180g" -R '+ref_g+'  -I '+bam+' --ERC GVCF --native-pair-hmm-threads '+threads+' --sample-ploidy 2 --min-pruning '+min_pruning+' --min-dangling-branch-length '+min_dangling+' -L '+CHR+' -O '+out_haplo+''
                    subprocess.Popen(haploCmd,shell=True).wait()

                else:
                    haploCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk HaplotypeCaller --java-options "-Xmx180g" -R '+ref_g+'  -I '+bam+' --ERC GVCF --native-pair-hmm-threads '+threads+' --sample-ploidy 2 --min-dangling-branch-length '+min_dangling+' -L '+CHR+' -O '+out_haplo+''
                    subprocess.Popen(haploCmd,shell=True).wait()

            else:

                if not (min_pruning == 'False'):
                    haploCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk HaplotypeCaller --java-options "-Xmx180g" -R '+ref_g+'  -I '+bam+' --ERC GVCF --native-pair-hmm-threads '+threads+' --sample-ploidy 2 --min-pruning '+min_pruning+' -L '+CHR+' -O '+out_haplo+''
                    subprocess.Popen(haploCmd,shell=True).wait()

                else:
                    haploCmd = 'module load tools java/1.8.0 gatk/4.1.8.1 && gatk HaplotypeCaller --java-options "-Xmx180g" -R '+ref_g+'  -I '+bam+' --ERC GVCF --native-pair-hmm-threads '+threads+' --sample-ploidy 2 -L '+CHR+' -O '+out_haplo+''
                    subprocess.Popen(haploCmd,shell=True).wait()


            # Generate sample map
            sample_map_file = out_dir+'/sample_map.'+CHR+'.txt'
            sample_map = open(sample_map_file,'a+')
            sample_map.write(bam_ID+'_'+CHR+'\t'+out_haplo+'\n')
            sample_map.close()
