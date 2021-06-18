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
    logi.write('The abundances of the non-MAG genes in the gene catalogue created by Prodigal 2.6.3, are obtained by mapping the reads\nnot included in the MAG set to the gene catalogue.\n\n')

# Inputs
# bam_files list
bam_files = glob.glob(bam_dir+'/*mapped.bam')

# annot files
annot_files = glob.glob(annot_dir+'/*-annotation.dmnd')
annot_files_str = ''
annot_IDs = list()
    # merge annotations
for annot_file in annot_files:
    annot_files_str += annot_file
    annot_IDs.append(annot_file.replace(annot_dir,'').replace('-annotation.dmnd',''))

if not glob.glob(annot_dir+'/*__annot.dmnd'):
    annot_db = annot_dir+'/'+'-'.join(annot_IDs)+'__annot.dmnd'

    if not (len(annot_IDs) == 1):
        mergeCmd='cat '+annot_files_str+' > '+annot_db+''    # merge the selected annotation dbs into one file
        subprocess.Popen(mergeCmd,shell=True).wait()
    else:
        mvCmd='mv '+annot_files_str+'  '+annot_db+''        # if only one annotation db, only rename
        subprocess.Popen(mvCmd,shell=True).wait()
else:
     annot_db = glob.glob(annot_dir+'/*__annot.dmnd')[0]

    # Create list of the genes that were successfully annotated by diamond
gene_annot__ids = {}
with open(annot_db,'r') as annot_data:
    for line in annot_data.readlines():
        (gene_ID,gene_annot,rest) = line.split('\t', 2) # keep two first fields of file
        gene_annot__ids[gene_ID.strip()] = gene_annot.strip()


# Will calculate total number of reads in each bam (mapped and unmapped)
# In case later user wants to get relative abundances
total_reads = out_dir+'/total_num_reads_BAMs.txt'
sample_list='Gene_Annot\tGene_ID\t'

# Index bam files
for bam in bam_files:
    if not os.path.isfile(bam+'.bai'):
        idxsamCmd='module load tools samtools/1.11 && samtools index '+bam+''
        #subprocess.Popen(idxsamCmd, shell=True).wait()

    sample = os.path.basename(bam).replace(ID+'.','').replace('.MAG_unmapped.bam','')
    all_genes_counts = out_dir+'/'+ID+'.'+sample+'.all_genes_counts.txt'

        #If the bam file has been indexed, continue
    if os.path.isfile(bam+'.bai'):
        if not os.path.exists(out_dir):
            mkdirCmd='mkdir -p '+out_dir+''
            subprocess.Popen(mkdirCmd,shell=True).wait()

        if not os.path.isfile(all_genes_counts):
                # extract total number of reads in bam file and append to common file
            totalCmd='module load tools samtools/1.11 && echo '+sample+' >> '+total_reads+' && samtools view -c '+bam+' >> '+total_reads+''
            subprocess.Popen(totalCmd,shell=True).wait()

                # calculate counts for all genes in .fna gene catalogue
            covCmd='module load tools samtools/1.11 && samtools idxstats '+bam+' | cut -f 1,3 > '+all_genes_counts+''
            subprocess.Popen(covCmd,shell=True).wait()


# Keep only genes successfully annotated by diamond from all genes
all_genes_files = glob.glob(out_dir+'/*all_genes_counts.txt')
annot_genes_files = list()

for file in all_genes_files:
    # file containing only annot
    sample = os.path.basename(file).replace(ID+'.','').replace('.all_genes_counts.txt','')
    sample_list += sample+'\t'
    annot_genes_counts = out_dir+'/'+ID+'.'+sample+'.annot_genes_counts.txt'
    annot_genes_files.append(annot_genes_counts)

    with open(file,'r') as all_genes_file, open(annot_genes_counts,'w+') as annot_genes:
        for line in all_genes_file.readlines():
            # if the given gene is found in the annot file keep it
            gene_ID = line.split('\t')[0].strip()
            if gene_ID in gene_annot__ids.keys():
                annot_genes.write(gene_annot__ids[gene_ID]+'\t'+line) # write the gene annotation + gene id + COUNTS
            else:
                pass


# 1 unique file per group with counts of annotates genes for all samples
all_counts_annot_genes = out_dir+'/'+ID+'.annot_counts_tmp.txt'
with open(all_counts_annot_genes,'w+') as final_annot_counts:
    final_annot_counts.write(sample_list+'\n')


pasteCmd='infiles="'+' '.join(annot_genes_files)+'" && cat '+annot_genes_files[0]+' | cut -f1,2 > GENEIDS && for i in $infiles; do sed -i -E "s/^.*\t.*\t//" $i; done && paste GENEIDS '+' '.join(annot_genes_files)+' >> '+all_counts_annot_genes+' && rm GENEIDS '+' '.join(annot_genes_files)+''
subprocess.Popen(pasteCmd,shell=True).wait()
# All annot genes files have the same genes, the total gene set. Thus, take first two columns (original gene ID, annotation) of the first file, and simply concatenate with all the
# counts in all files.
