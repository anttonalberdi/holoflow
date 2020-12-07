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
  - ***metagenomics_DR.py***  - Dereplication of metagenomic bins produced by either *metagenomics_IB* or *metagenomics_CB*. 
  
  
  
These are designed to be called from the command line and require the following arguments (**{only in PREPROCESSING}**,**[optional arguments]**):  
```bash
  -f INPUT            File containing input information.
  -d WORK_DIR         Output directory.
  -t THREADS          Thread maximum number to be used by Snakemake.  
  {-r REF_GENOME}     Reference genome(s) file path to be used in read mapping.
  [-k KEEP_TMP]       If present, keep temporal directories - NOT IN PREPAREGENOMES.
  [-l LOG]            Desired pipeline log file path.
  [-c CONFIG]         Configuration file full path.
  
```  

 
#### Config files description
A template *config.yaml* file can be found in every workflow directory. 

#### Input files description
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
  2. Coassembly group.  
  3. Original full path/name of **FORWARD** input file.  
  4. Original full path/name of **REVERSE** input file.  
   * Optimally the metagenomic .fastq files would come from PPR_03-MappedToReference, the last preprocessing step.
  
- Example:

|   |   |   |   |
| --- | --- | --- | --- |
| Sample1 | CoassemblyGroup1 | /home/Sample1_1.fq | /home/Sample1_2.fq |  
| Sample2 | CoassemblyGroup2 | /home/Sample2_1.fq | /home/Sample1_2.fq |  
| Samplen | CoassemblyGroup3 | /home/Samplen_1.fq | /home/Samplen_2.fq |
  

##### *metagenomics_DR.py*

  1. Coassembly group or sample group name.  
  2. Input directory path where all *.fa* bins to dereplicate are.
  
- Example:

|   |   |   |
| --- | --- | --- |
| GroupA | /home/directory_samplesA |
| GroupB | /home/directory_samplesB |

  
 
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

- Config file *config.yaml*, in which the user may be interested to customise:
  1. Quality filtering - specific adapter sequences, minimum quality, character separating the mate read number.


#### Metagenomics - Individual Assembly & Coassembly
- *Snakefile* - which contains rules for:
  1. Metagenomic assembly using **metaSpades** or **megahit**
  2. Read mapping to assembly using **bwa mem** 
  3. Contig binning using **Metabat**, **MaxBin** (and **Concoct** #### NOT YET)
  4. Binner result integration using **DasTool** 
  
- Config file *config.yaml*, in which the user may be interested to customise:
  1. Assembler - choose between the mentioned options by writing *megahit* or *spades*
  2. Minimum contig length - minimum bp per contig in final assembly file.

  
#### Metagenomics - Dereplication
- *Snakefile* - which contains rules for:
  1. Bin Dereplication using **dRep**
  2. Bin assembly improvement (contig elongation and scaffolding) using SSPACE. ##### UNDER CONSTRUCTION
  3. Phylogenetic analysis and taxonomic assignation ##### UNDER CONSTRUCTION 
  
- Config file *config.yaml*, in which the user may be interested to customise:
  1. Desired contig scaffolding or not, by setting SSPACE *True/False*



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
python ${holoflowpath}/preprocessing.py -f ${projectpath}/input.txt -d ${projectpath}/workdir -r ${projectpath}/reference_genomes.fna -c ${projectpath}/config.yaml -l ${projectpath}/log_file.log -t 40
```

- *job execution* in Computerome2 example:
```bash
 qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${projectpath}/job_error_file.err -o ${projectpath}/job_out_file.out -l nodes=1:ppn=40,mem=180gb,walltime=5:00:00:00 -N JOB_ID ${projectpath}/first_job_preprocessing.sh

```
  Note that the job parameters: *ppn*, *nodes*, *memory*, *wall time* ... can and ought to be customised optimally for every job type.





