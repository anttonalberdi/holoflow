# holoflow
Bioinformatics pipeline for hologenomics data generation and analysis

Snakemake is a workflow management system which requires from a *Snakefile* and a *config* file. This is a Bioinformatics pipeline implemented with Snakemake.

## Files and directories
### Main directory

The main *holoflow* directory contains a given number of Python scripts which work as launchers for the different **workflow programs** in the pipeline:

  - ***preparegenomes.py***   - Merge all potential reference genomes to sample into a single *.fna* file to be used in preprocessing.py.  
  - ***preprocessing.py***    - Data preprocessing from quality to duplicate sequences for further downstream analysis.
  - ***metagenomics_IB.py***  - Individual assembly-based analysis and metagenomics binning. 
  - ***metagenomics_CB.py***  - Coassembly-based analysis and metagenomics binning. 
  - ***metagenomics_DR.py***  - Dereplication and Annotation of metagenomic bins produced by either *metagenomics_IB* or *metagenomics_CB*. 
  - ***metagenomics_FS.py***  - Final statistical report of dereplicated bins obtained with *metagenomics_DR.py*. 
  - ***metagenomics_AB.py***  - Functional annotation of (co-)assembly file with DRAM.  
  - ***genomics.py***         - Variant calling, Phasing (for HD) and Imputation (for LD) with *genomics.py*. 

  
  
These are designed to be called from the command line and require the following arguments:  
```bash
REQUIRED ARGUMENTS:
  -f INPUT            File containing input information.
  -d WORK_DIR         Output directory.
  -t THREADS          Thread maximum number to be used by Snakemake.
  -W REWRITE          Wants to re-run the worfklow from scratch: remove all directories previous runs. - NOT IN PREPAREGENOMES.
  -g REF_GENOME   Reference genome(s) file path to be used in read mapping. Unzipped for genomics. - only in PREPROCESSING, GENOMICS.  
  -adapter1 ADAPTER1 Adapter sequence 1 for removal. - only in PREPROCESSING.
  -adapter2 ADAPTER2 Adapter sequence 2 for removal. - only in PREPROCESSING. 
  -Q DATA QUALITY]     Low depth (LD) or High depth (HD) data set. - only in GENOMICS.
  -vc VAR CALLER       Variant caller to choose: 1 {bcftools/samtools}, 2 {GATK}, 3 {ANGSD}. - only in GENOMICS.
  -N JOB ID            ID of the sent job, so another different-N-job can be run simultaneously. - only in GENOMICS, METAGENOMICS IB, AB.

OPTIONAL ARGUMENTS:
  -r REF_PANEL        Reference panel necessary for likelihoods update and imputation of LD variants. - only in GENOMICS.
  -k KEEP_TMP         If present, keep temporal directories - NOT IN PREPAREGENOMES.
  -l LOG              Desired pipeline log file path.
  -c CONFIG           Configuration file full path.
  
```  
 
 
### Config files description
A template *config.yaml* file can be found in every workflow directory. 

### Input files description
A template *input.txt* file can be found in every workflow directory.  
See *input.txt* file description for every workflow:
In all cases, columns must be delimited by a simple space and **no blank lines should be found in the end of the file**.  
Those lines starting by # won't be considered.  
  
##### *preparegenomes.py*

  1. Reference genomes ID. **No spaces or undersquares** between different words in identifier.  
  2. Reference genome full path/name.   
  3. Desired output data base with all genomes name. **No spaces**, undersquares or other separators allowed. *All those reference genomes which should be in the same DB should have the same ID in this field*.  
  
  **The fields 1 and 3 must be different**
  
- Example:  

*Heads-up*: you can generate more than one DB at a time for different projects, be aware that preprocessing only takes ONE DB at a time with all reference genomes to be mapped to a set of samples in a given project.

|   |   |   |
| --- | --- | --- |
| Genomeone   | /home/Genomeone.fq      | DBone  |
| Genometwo   | /home/Genometwo.fq.gz   | DBtwo  |
| Genomethree | /home/Genomethree.fq    | DBone  |
| Genomen     | /home/Genomen.fq        | DBn    |


##### *preprocessing.py*  &  *metagenomics_IB.py*

  1. Sample name.  
  2. Original full path/name of **FORWARD** input file. This can be both *.gz* or not compressed.  
  3. Original full path/name of **REVERSE** input file. This can be both *.gz* or not compressed.  
  
- Example:

|   |   |   |   |
| --- | --- | --- | --- |
| Sample1 | /home/Sample1_1.fq | /home/Sample1_2.fq |  
| Sample2 | /home/Sample2_1.fq | /home/Sample1_2.fq |  
| Samplen | /home/Samplen_1.fq | /home/Samplen_2.fq |  


##### *metagenomics_CB.py*

  1. Sample name.  
  2. Coassembly group: **assumed to be the same as in preprocessing -N job if preprocessing has been run (PPR_03-MappedToReference job directory ID)**.  
  3. Original full path/name of **FORWARD** input file.  
  4. Original full path/name of **REVERSE** input file.  
Optimally the metagenomic .fastq files would come from PPR_03-MappedToReference, the last preprocessing step.
  
- Example:

|   |   |   |   |
| --- | --- | --- | --- |
| Sample1 | CoassemblyGroup1 | /home/Sample1_1.fq | /home/Sample1_2.fq |  
| Sample2 | CoassemblyGroup2 | /home/Sample2_1.fq | /home/Sample1_2.fq |  
| Samplen | CoassemblyGroup3 | /home/Samplen_1.fq | /home/Samplen_2.fq |
  

##### *metagenomics_DR.py*

  1. Coassembly group or sample group name.  
  2. Input directory path where all *.fa* bins to dereplicate and the respective *ID*_DASTool_summary.txt files are.
  
- Example:

|   |   |   |
| --- | --- | --- |
| GroupA | /home/directory_samplesA |
| GroupB | /home/directory_samplesB |
  

##### *metagenomics_FS.py*

  1. Coassembly group or sample group name.  
  2. Input directory path where the group's/samples' in the group original metagenomic *_1.fastq* & *_2.fastq* files are.
  3. Input directory path where all dereplicated *.fa* bins are.
  4. Input directory path where .gff annotation files respective to each dereplicated bin is found.
  
- Example:

|   |   |   |   |
| --- | --- | --- | --- |
| DrepGroup1 | /home/PPR_03-MappedToReference/DrepGroup1 | /home/MDR_01-BinDereplication/DrepGroup1/dereplicated_genomes | /home/MDR_02-BinAnnotation/DrepGroup1/bin_funct_annotations |
| DrepGroup2 | /home/PPR_03-MappedToReference/Sample1 | /home/MDR_01-BinDereplication/Sample1/dereplicated_genomes | /home/MDR_02-BinAnnotation/DrepGroup2/bin_funct_annotations | 
| DrepGroup2 | /home/PPR_03-MappedToReference/Sample2 | /home/MDR_01-BinDereplication/Sample2/dereplicated_genomes | /home/MDR_02-BinAnnotation/DrepGroup3/bin_funct_annotations |


##### *metagenomics_AB.py*

  1. (Co-)Assembly or group ID.   
  2. Path to assembly file.  
  
- Example:

|   |   |   |
| --- | --- | --- |
| GroupA | /home/dir/assembly_A.fa |
| GroupB | /home/second/dir/assembly_B.fna.gz |


##### *genomics.py*

  1. Sample group name to analyse.  
  2. Path to directory containing host reads BAM alignment sorted files - If *preprocessing.py* was used, these are the resulting *ref* BAMs path.   
  3. Chromosome list. This should be a text file with a single column depicting chromosome IDs. Note that **the given chromosome IDs should be in accordance with the provided reference genome**, otherwise these won't be detected by Holoflow.  
  
- Example:  

|   |   |   |
| --- | --- | --- |
| Chicken_samples   | /home/path/to/chicken/bams      |  /home/path/to/chicken_chrlist.txt  |
| Cervid_samples   | /home/path/to/cervid/PPR_03-MappedToReference   | /home/path/to/cervid_chrlist.txt  |
| Cavia_samples | /home/path/to/cavia/bams     | /home/path/to/cavia_chrlist.txt  |


 
### Workflows - Specific directories

#### Preparegenomes
- *Snakefile* - Continuing *preparegenomes.py*'s job, which takes as input the full paths of the given reference genomes, reformats its read IDs and merges them into a single *data_base.fna* file, the *Snakefile* contains rules for:  
  1. Indexing the resulting DB using **bwa** and **samtools**
  2. Compressing the full set of DB-related files into a *data_base.tar.gz* file.


#### Preprocessing
- *Snakefile* - which contains rules for:
  1. Quality filtering using **AdapterRemoval**
  2. Duplicate read removal using **seqkit rmdup**
  3. Mapping reads against reference genome(s) using **bwa mem**

- Config file *config.yaml*, in which the user may be interested in customising:
  1. Quality filtering - specific adapter sequences, minimum quality, character separating the mate read number.


#### Metagenomics - Individual Assembly & Coassembly
- *Snakefile* - which contains rules for:
  1. Metagenomic assembly using **megahit**. In Individual Assembly also **metaSpades** available.  
  2. Read mapping to assembly using **bwa mem** 
  3. Contig binning using **Metabat**, **MaxBin**. In Coassembly also binning by **Concoct**.  
  4. Binner result integration using **DasTool** 
  
- Config file *config.yaml*, in which the user may be interested in customising:
  1. Assembler - choose between the mentioned options by writing *megahit* or *spades*
  2. Minimum contig length - minimum bp per contig in final assembly file.

  
#### Metagenomics - Dereplication
- *Snakefile* - which contains rules for:
  1. Bin Dereplication using **dRep**.
  2. Bin Gene Annotation with **Prokka**.
  3. Bin Taxonomic Classification with **GTDB-Tk**.
  4. Obtain GTDB phylogenetic subtree of MAGs.

  
#### Metagenomics - Final Statistics
- *Snakefile* - which contains rules for:
  1. Mapping metagenomic reads to dereplicated MAGs - number and % of mapped reads.
  2. Obtaining coverage statistics of contigs and MAGs in used samples.
  3. Retrieve quality statistics (CheckM) and summary plot of the MAGs.
  4. Get coverage of KEGG KO single-copy core genes in MAGs. 

#### Metagenomics - Assembly Based
- *Snakefile* - which contains rules for:
  1. DRAM functional annotation and distilling of an assembly file.   
  
#### Genomics
- *Snakefile* - which contains rules for:  
 a. Variant calling with **BCFtools**, **GATK** or **ANGSD** (## Latter UNDER CONSTRUCTION ##)  

  -> *High depth samples*  
 b. Filtering with **BCFtools** or **GATK**  
 c. Phasing with **shapeit4**  

  -> *Low depth samples*  
 b. Likelihoods update with **Beagle** using a high-depth reference panel  
 c. Genotype imputation with **Beagle**   
  
- Config file *config.yaml*, in which the user may be interested in customising:
  1. Choose between HD - for high depth seqs OR LD - for low depth seqs.
  2. Variant calling - BCFtools
    - mpileup
      * Coefficient for downgrading mapping quality for reads containing excessive mismatches - *degr_mapp_qual*. Default 50.  
      * Minimum mapping quality - *min_mapp_qual*. Default to 0.  
      * Minimum base quality - *min_base_qual*. Default to 13.  
      * Specific chromosome region. Default False.  
    - call
      * Multicaller mode: alternative model for multiallelic and rare-variant calling designed to overcome known limitations. 
      * Keep only variants and not indels. 
      
  3. Variant calling - GATK 
      * Parameters to obtain more agressive variants: *min_pruning* and *min_dangling*.
   
  4. Variant calling - ANGSD
      * Choose model (1/2) between samtools or GATK.
      * Output log genotype likelihoods to a file or not.
      * How to estimate minor and major alleles (1/2): 1 = from likelihood data ; 2 = from count data.
      * Estimate posterior genotype probability based on the allele frequency as a prior (True/False).
  5. HD Filtering - BCFtools 
      * Quality of SNPs that want to be kept. Default to 30.
  6. HD Filtering - GATK
      * Quality of SNPs that want to be kept. Default to 30.
      * QD: Quality by depth. Find more information [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).
      * FS: Fisher strand. Find more information [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).

  7. HD Phasing 
      * --geno filters out all variants with missing call rates exceeding the provided value to be removed. Default to 0.
      * Provide a Genetic map. Default to False, else provide path.


## Usage in Computerome

### Get started: download Holoflow repository
Clone the repository by running the following command on your command line:

```bash
git clone -b nurher --single-branch https://github.com/anttonalberdi/holoflow.git
```

### Execute Holoflow *.py* workflow launchers
These should be **executed as jobs**, therefore a *.sh* script should be generated which will call the desired Holoflow workflow:

- *.sh* example script for *preprocessing.py* called ***first_job_preprocessing.sh***:
```bash
#Declare full path to the project directory (the .sh file will be stored here as well)
projectpath=/full/path/project1
#Declare full path to holoflow
holoflowpath=/full/path/holoflow
#Run holoflow
python ${holoflowpath}/preprocessing.py -f ${projectpath}/input.txt -d ${projectpath}/workdir -g ${projectpath}/reference_genomes.fna -adapter1 'ATGCT' -adapter2 'CTTGATG' -c ${projectpath}/config.yaml -l ${projectpath}/log_file.log -t 40 -N First_job
```

- *job execution* in Computerome2 example:
```bash
 qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${projectpath}/job_error_file.err -o ${projectpath}/job_out_file.out -l nodes=1:ppn=40,mem=180gb,walltime=5:00:00:00 -N JOB_ID ${projectpath}/first_job_preprocessing.sh

```
  Note that the job parameters: *ppn*, *nodes*, *memory*, *wall time* ... can and ought to be customised optimally for every job type.





