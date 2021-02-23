#22.11.2020 - Holoflow 0.1.

import subprocess
import argparse
import os
import sys
import glob
import time
import gzip
import numpy as np


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
annot_dir=args.annot_dir
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

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

    # Prepare mag, bam data and ID
    mag_list=glob.glob(str(mag_dir)+'/*.fa')
    bam_list=glob.glob(str(bam_dir)+'/*.bam')
    gff_list = glob.glob(annot_dir+'/*.gff')

    for i in range(len(mag_list)):
        mag = mag_list[i]
        mag_ID = os.path.basename(mag).replace('.fa','')

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


        for bam in bam_list:
            sample = os.path.basename(bam).replace('.bam','')
            new_bam = out_dir+'/'+mag_ID+'_'+sample+'.bam'
            sample_counts_tmp = out_dir+'/'+mag_ID+'_'+sample+'.counts.txt'

            if not os.path.isfile(new_bam):
            # Split bams into MAGs
            # Now BAM headers are only the contig ID - Removed MAG_ID-
                samtoolsCmd='module load tools samtools/1.11 && samtools view -h '+bam+' | grep "'+mag_ID+'-" | sed "s/'+mag_ID+'-//" | samtools view -bS - | htseq-count -t CDS -r pos -f bam - '+gtf+' > '+sample_counts_tmp+''
                subprocess.Popen(samtoolsCmd,shell=True).wait()

            else:
                htseqCountsCmd='module load tools && htseq-count -t CDS -r pos -f bam '+new_bam+' '+gtf+' > '+sample_counts_tmp+'' ## ?? --nonunique all ??
                subprocess.Popen(htseqCountsCmd,shell=True).wait()


    #Some files will be empty -> remove them
    try:
        rmCmd='find '+out_dir+' -size 0 -delete'
        subprocess.Popen(rmCmd,shell=True).wait()
    except:
        pass

    ## Handle coverage and IDs

    # Read KO_db into a dictionary [Uniprot]=KO
    with gzip.open(KO_db,'rt') as kos_db:
        KO_database = {}
        for line in kos_db:
            (key,val) = line.split()
            KO_database[key] = val


    ## Get coverage of annotated genes
    for mag in mag_list:
        sample_list = 'KO\t'
        KO_times = {}
        n = 0

        mag_ID = os.path.basename(mag).replace('.fa','')
        mag_annot = annot_dir+'/'+mag_ID+'.gtf'
        mag_counts_tmp = out_dir+'/'+mag_ID+'_counts_tmp.txt'

        counts_list = glob.glob(out_dir+'/'+mag_ID+'_*.counts.txt')
        counts_string = ''
        for file in counts_list:
            counts_string+=file.strip()+' '
            sample = os.path.basename(file).replace('.counts.txt','').replace(mag_ID+'_','')
            sample_list+=sample+'\t'

        pasteCmd='infiles="'+counts_string+'" && for i in $infiles; do sed -i -E "s/^.*\t//" $i; done && cut -f1 '+counts_list[0]+' > UNIPROT && paste UNIPROT '+counts_string+' > '+mag_counts_tmp+' && rm UNIPROT'
        subprocess.Popen(pasteCmd,shell=True).wait()



        mag_counts = out_dir+'/'+mag_ID+'_counts.txt'
    # Reformat - Translate annotation in counts file UniProt -> KO
        with open(mag_counts_tmp,'r') as tmp_counts, open(mag_counts,'w+') as final_counts:
            final_counts.write(sample_list+'\n')

            for line in tmp_counts.readlines():
                line=line.split('\t',1) # max number of splits 1
                uniprot=line[0]
                counts=line[1]

                try:
                    KO = KO_database[str(uniprot).strip()]
                    # Write new data to final counts
                    final_counts.write(KO+'\t'+counts)

                    ## Generate file ONLY for KO counts in the list
                    with open(KO_genes,'r') as ko_genes:
                        for line in ko_genes.readlines():
                            if KO in line:
                            # Write new data to ko counts
                                if not KO in KO_times.keys():
                                    KO_times[KO] = []
                                    KO_times[KO].append(counts.split('\t'))
                                else:
                                    KO_times[KO].append(counts.split('\t'))
                except:
                    pass


        KO_counts = out_dir+'/'+mag_ID+'_KO_counts.txt'
        with open(KO_counts,'w+') as ko_counts:
            sample_list = sample_list.split('\t')[:-1]
            sample_list.insert(len(sample_list),'N')
            sample_list = ('\t').join(sample_list)
            ko_counts.write(sample_list+'\n')

            for key in KO_times.keys():
                n = len(KO_times[key])
                counts_sum = np.array(KO_times[key]).astype(int)
                counts_sum = np.sum(counts_sum,axis=0)
                counts_sum = counts_sum.tolist()
                counts_sum = '\t'.join(str(v) for v in counts_sum)

                ko_counts.write(key+'\t'+str(counts_sum)+'\t'+str(n)+'\n')



    #os.remove(mag_counts_tmp)
