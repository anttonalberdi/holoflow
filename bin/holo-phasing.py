## 26.02.21 - Holoflow 0.1
import subprocess
import argparse
import os
import glob
import time


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-filt_dir', help="filtered variants directory", dest="filt_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-chr_list', help="chromosome list file path", dest="chr_list", required=True)
parser.add_argument('-geno', help="number of missing genotypes allowed", dest="geno", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-gmap', help="gmap", dest="gmap", required=True)
args = parser.parse_args()


filt_dir=args.filt_dir
out_dir=args.out_dir
chr_list=args.chr_list
geno=args.geno
ID=args.ID
log=args.log
threads=args.threads
gmap=args.gmap


## Run
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Write to log
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    with open(str(log),'a+') as logi:
        logi.write('\t\t'+current_time+'\tPhasing of HD data - '+ID+'\n')
        logi.write(' \n\n')

    chromosome_list = list()
    with open(chr_list,'r+') as chr_data:
        for chr in chr_data.readlines():
            chromosome_list.append(chr.strip())

    for CHR in chromosome_list:
        input = filt_dir+'/'+ID+'.HD_SNPs_'+CHR+'.vcf.gz'
        plink_tmp_output_base = out_dir+'/'+ID+'.plink_tmp.HD_SNPs_'+CHR
        plink_output_base = out_dir+'/'+ID+'.plink.HD_SNPs_'+CHR
        output = out_dir+'/'+ID+'_'+CHR+'.filt_phased.vcf.gz'

        # Plink filtration of SNPs before phasing
        plink1Cmd='module load plink2/1.90beta6.17 && plink --vcf '+input+' --double-id --make-bed --allow-extra-chr --keep-allele-order  --real-ref-alleles --set-missing-var-ids "@:#\$1,\$2" --out '+plink_tmp_output_base+''
        subprocess.Popen(plink1Cmd,shell=True).wait()

        plink2Cmd='module load plink2/1.90beta6.17 && plink --bfile '+plink_tmp_output_base+' --double-id --allow-extra-chr --keep-allele-order  --real-ref-alleles --geno '+geno+' --recode vcf-iid bgz --out '+plink_output_base+''
        subprocess.Popen(plink2Cmd,shell=True).wait()

    # Filter output
        if not os.path.isfile(plink_output_base+'.vcf.csi'):
            indexCmd='module load bcftools/1.11 && bcftools index --threads '+threads+' '+plink_output_base+'.vcf.gz'
            subprocess.Popen(indexCmd,shell=True).wait()


        if not (gmap == 'False'):
            phasingCmd= 'module load shapeit4/4.1.3 && shapeit4 --input '+plink_output_base+'.vcf.gz --map '+gmap+' --region '+CHR+' --thread '+threads+' --output '+output+' --sequencing'
            subprocess.Popen(phasingCmd,shell=True).wait()

        else:
            phasingCmd= 'module load shapeit4/4.1.3 && shapeit4 --input '+plink_output_base+'.vcf.gz --region '+CHR+' --thread '+threads+' --output '+output+' --sequencing'
            subprocess.Popen(phasingCmd,shell=True).wait()

        # Index phased panel
        idxCmd='module load tabix/1.2.1 && tabix '+output+''
        subprocess.Popen(idxCmd,shell=True).wait()


    # Concatenate all CHR phased files into one ref panel
    ref_panel_phased = out_dir+'/'+ID+'_RefPanel-Phased.vcf.gz'
    phased_files = glob.glob(out_dir+'/'+ID+'_*filt_phased.vcf.gz')
    files_to_concat = out_dir+'/'+ID+'_files_to_concat.txt'
    with open(files_to_concat,'w+') as concat:
        for file in phased_files:
            concat.write(file.strip()+'\n')

    # make sure chr in same order chr list
    concatCmd= 'module load bcftools/1.11  && bcftools concat -f '+files_to_concat+' -Oz -o '+ref_panel_phased+' && rm '+files_to_concat+''
    subprocess.Popen(concatCmd,shell=True).wait()
