#22.11.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import sys
import glob
import time
import gzip


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-bam_dir', help="input bam from mapped MAGs to .fastq directory", dest="bam_dir", required=True)
parser.add_argument('-mag_dir', help="originally dereplicated mags", dest="mag_dir", required=True)
parser.add_argument('-annot_dir', help="annotation directory", dest="annot_dir", required=True)
parser.add_argument('-out_dir', help="main output directory", dest="out_dir", required=True)
parser.add_argument('-KO_db', help="data base UniProt-KO", dest="KO_db", required=True)
parser.add_argument('-KO_list', help="KO genes to find", dest="KO_genes", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

bam_dir=args.bam_dir
mag_dir=args.mag_dir
annot_dir=srgs.annot_dir
out_dir=args.out_dir
KO_db=args.KO_db
KO_genes=args.KO_genes
ID=args.ID
log=args.log
threads=args.threads



# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\t\t'+current_time+'\t - '+ID+'\n')
    logi.write('\t')


# Prepare mag, bam data and ID
mag_list=glob.glob(str(mag_dir)+'/*.fa')
bam_list=glob.glob(str(bam_dir)+'/*.bam')
gff_list = glob.glob(annot_dir+'/*.gff')

for i in range(len(mag_list)):
    mag = mag_list[i]
    mag_ID = os.path.basename(mag).replace('.fa','')


    for bam in bam_list:
        sample = os.path.basename(bam).replace('.bam','')
        new_bam = out_dir+'/'+mag_ID+'_'+sample+'.bam'

        if not os.path.isfile(new_bam):
        # Split bams into MAGs
        # Now BAM headers are only the contig ID - Removed MAG_ID-
            samtoolsCmd='module load tools samtools/1.11 && samtools view -h '+bam+' | grep "'+mag_ID+'-" | sed "s/'+mag_ID+'-//" | samtools view -bS - > '+new_bam+''
            subprocess.Popen(samtoolsCmd,shell=True).wait()

    # Reformat GFF > GTF
    gff = gff_list[i]
    gtf = gff.replace('.gff','.gtf')
    tmp_prokka = gff.replace('.gff','_tmp_prokka')
    tmp_uniprot = gff.replace('.gff','_tmp_uniprot')


        # retrieve current directory
    file = os.path.dirname(sys.argv[0])
    curr_dir = os.path.abspath(file)

    gtfCmd='bash '+curr_dir+'/holo-create_gtf.sh '+gff+' > '+gtf+''
    subprocess.Popen(gtfCmd,shell=True).wait()


# Some bam files will be empty -> remove them
try:
    rmCmd='find '+out_dir+' -size 0 -delete'
    subprocess.Popen(rmCmd,shell=True).wait()
except:
    pass



## Handle coverage and IDs

# Read KO_db into a dictionary [Uniprot]=KO
with gzip.open(KO_db,'r') as kos_db:
    KO_database = {}
    for line in kos_db:
        (key,val) = line.split()
        KO_database[key] = val


sample_list = 'KO '
## Get coverage of annotated genes
for mag in mag_list:
    mag_ID = os.path.basename(mag).replace('.fa','')
    mag_annot = annot_dir+'/'+mag_ID+'.gtf'
    mag_counts_tmp = out_dir+'/'+mag_ID+'_counts.txt_tmp'

    mag_bams_list = glob.glob(out_dir+'/'+mag_ID+'_*.bam')
    mag_bams = ''
    for bam in mag_bams_list:
        mag_bams+=bam+' '
        sample = os.path.basename(bam).replace('.bam','')
        sample_list+=sample+' '

    htseqCountsCmd='module load tools && htseq-count -t CDS -r pos -f bam '+mag_bams+' '+mag_annot+' > '+mag_counts_tmp+'' ## ?? --nonunique all ??
    subprocess.Popen(htseqCountsCmd,shell=True).wait()


## Reformat - Translate annotation in counts file UniProt -> KO
    mag_counts = out_dir+'/'+mag_ID+'_counts.txt'
    KO_counts = out_dir+'/'+mag_ID+'_KO_counts.txt'
    with open(mag_counts_tmp,'r+') as tmp_counts, open(mag_counts,'w+') as final_counts, open(KO_counts,'w+') as ko_counts:
        data = tmp_counts.readlines()
        final_counts.write(sample_list+'\n')

        for line in data:
            line=line.split('\t',1) # max number of splits 1
            uniprot=line[0]
            counts=line[1]
            KO = KO_database[str(uniprot).strip()]
            print(KO)

            # Write new data to final counts
            final_counts.write(KO+'\t'+counts+'\n')


## Generate file ONLY for KO counts in the list
            with open(KO_genes,'r') as ko_genes:
                ko_counts.write(sample_list+'\n')
                if str(KO).strip() in ko_genes:
                    # Write new data to ko counts
                    ko_counts.write(KO+'\t'+counts+'\n')
