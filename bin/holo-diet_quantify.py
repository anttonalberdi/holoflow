#06.05.2021 - Holoflow 0.1.
import subprocess
import argparse
import os
import time
import glob


#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-annot_dir', help="annotation directory ", dest="annot_dir", required=True)
parser.add_argument('-bam_dir', help="directory with mappings .fq + GC", dest="bam_dir", required=True)
parser.add_argument('-out_dir', help="out_dir", dest="out_dir", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-ID', help="ID", dest="ID", required=True)
parser.add_argument('-log', help="pipeline log file", dest="log", required=True)
args = parser.parse_args()


annot_dir=args.annot_dir
bam_dir=args.bam_dir
out_dir=args.out_dir
t=args.threads
ID=args.ID
log=args.log



# Run
# Write to log
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
with open(str(log),'a+') as logi:
    logi.write('\tHOLOFLOW\tMETAGENOMICS\n\t\t'+current_time+'\t - '+ID+'\n')
    logi.write('   \n\n')

# Inputs
# annot file
annot_file = glob.glob(annot_dir+'/*-annotation.dmnd')[0]
# bam_files list
bam_files = glob.glob(bam_dir+'/*mapped.bam')


    # Create list of the genes that were successfully annotated by diamond
gene_annot__ids = {}
with open(annot_file,'r') as annot_data:
    for line in annot_data.readlines():
        (gene_ID,gene_annot) = line.split('\t', 1) # keep two first fields of file
        gene_annot__ids[gene_ID.strip()] = gene_annot.strip()


# Will calculate total number of reads in each bam (mapped and unmapped)
# In case later user wants to get relative abundances
total_reads = out_dir+'/total_num_reads_BAMs.txt'
sample_list='Gene_ID\t'

# Index bam files
for bam in bam_files:
    if not os.path.isfile(bam+'.bai'):
        idxsamCmd='module load tools samtools/1.11 && samtools index '+bam+''
        subprocess.Popen(idxsamCmd, shell=True).wait()

    sample = os.path.basename(bam).replace('bam_dir','').replace('.mapped.bam','')
    sample_list += sample+'\t'
    all_genes_counts = out_dir+'/'+ID+'.'+sample+'.all_genes_counts.txt'

        # If the bam file has been indexed, continue
    if os.path.isfile(bam+'.bai'):
        if not os.path.isfile(all_genes_counts):
                # extract total number of reads in bam file and append to common file
            totalCmd='module load tools samtools/1.11 && echo '+sample+' >> '+total_reads+' && samtools view -c '+bam+' >> '+total_reads+''
            subprocess.Popen(totalCmd,shell=True).wait()

                # calculate counts for all genes in .fna gene catalogue
            covCmd='module load tools samtools/1.11 && samtools idxstats '+bam+' | cut -f 1,3 > '+all_genes_counts+''
            subprocess.Popen(covCmd,shell=True).wait()


# Keep only genes successfully annotated by diamond from all genes
all_genes_files = glob.glob(out_dir+'/*all_genes_counts.txt')

for file in all_genes_files:
        # file containing only annot
     annot_genes_counts = out_dir+'/'+ID+'.'+sample+'.annot_genes_counts.txt'

     with open(file,'r') as all_genes_file, open(annot_genes_counts,'w+') as annot_genes:
         for line in all_genes_file.readlines():
             # if the given gene is found in the annot file keep it
             gene_ID = line.split()[0].strip()
             if gene_ID in gene_annot__ids.keys():
                 annot_genes.write(gene_annot__ids[gene_ID]+'\t'+line) # write the gene annotation + gene id + COUNTS
             else:
                 pass


# Merge counts of all samples in one file
annot_genes_files = glob.glob(out_dir+'/*all_genes_counts.txt')
annot_genes_files_string = ''
for file in annot_genes_files:
    annot_genes_files_string += file+' '

# 1 unique file per group with counts of annotates genes for all samples
all_counts_annot_genes = out_dir+'/'+ID+'.annot_counts_tmp.txt'

pasteCmd='infiles="'+annot_genes_files_string+'" && cat '+annot_genes_files[0]+' | cut -f1,2 > UNIPROT && for i in $infiles; do sed -i -E "s/^.*\t.*\t//" $i; done && paste UNIPROT '+annot_genes_files_string+' > '+all_counts_annot_genes+' && rm UNIPROT'
subprocess.Popen(pasteCmd,shell=True).wait()
