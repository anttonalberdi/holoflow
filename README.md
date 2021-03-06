# holoflow
Bioinformatics pipeline for hologenomics data generation and analysis

Snakemake is a workflow management system which requires from a *Snakefile* and a *config* file. This is a Bioinformatics pipeline for hologenomics data generation and analysis implemented with Snakemake.

## Files and directories
### Main directory

The main *holoflow* directory contains a given number of Python scripts which work as launchers for the different **workflow programs** in the pipeline:

  - *preparegenomes.py*   - Merge all potential reference genomes to sample into a single *.fna* file to be used in preprocessing.py.  
  - *preprocessing.py*    - Data preprocessing from quality to duplicate sequences for further downstream analysis.
  - *metagenomics_IA.py*  - Individual assembly-based assembly and metagenomics binning. 
  
  
These are designed to be called from the command line and require the following arguments ([optional arguments]):  
```bash
  -f INPUT            File containing input information.
  -d WORK_DIR         Output directory.
  -t THREADS          Thread maximum number to be used by Snakemake.  
  [-l LOG]            Desired pipeline log file path.
  [-c CONFIG]         Configuration file full path.
  
```  
  
#### Input files description
Find *input.txt* file description for every workflow.  
In all cases, columns must be delimited by a simple space and no blank lines should be found in the end of the file.  
Those lines starting by # won't be considered.  
  
##### *preparegenomes.py*

  1. Reference genomes ID. **No spaces or undersquares** between different words in identifier.  
  2. Reference genome full path/name.   
  3. Desired output data base with all genomes name. **No spaces**, undersquares or other separators allowed. *All those reference genomes which should be in the same DB should have the same ID in this field*.  
  
- Example:

|   |   |   |
| --- | --- | --- |
| Genomeone   | /home/Genomeone.fq      | DBone  |
| Genometwo   | /home/Genometwo.fq.gz   | DBtwo  |
| Genomethree | /home/Genomethree.fq    | DBone  |
| Genomen     | /home/Genomen.fq        | DBn    |


##### *preprocessing.py*  &  *metagenomics_IA.py*

  1. Sample name.  
  2. Assembly group (If not *metagenomics/coassembly* this field will be ignored - nevertheless, it is important that is not omitted when writing the input file).   
  3. Original full path/name of input file/s. These can be both *.gz* or not compressed.  
  
- Example:

|   |   |   |
| --- | --- | --- |
| Sample1 | Group1 | /home/Sample1_1.fq |
| Sample1 | Group1 | /home/Sample1_2.fq |
| Sample2 | Group1 | /home/Sample2_1.fq |
| Sample2 | Group1 | /home/Sample1_2.fq |
| Samplen | Groupn | /home/Samplen_1.fq |
| Samplen | Groupn | /home/Samplen_2.fq |
  
  
  
 
### Workflows - Specific directories

#### Preparegenomes
- *Snakefile* - Continuing *preparegenomes.py*'s job, which takes as input the full paths of the given reference genomes, reformats its read IDs and merges them into a single *data_base.fna* file, the *Snakefile* contains rules for:  
  1. Indexing the resulting DB using **bwa** and **samtools**
  2. Compressing the full set of DB-related files into a *data_base.fna.tar.gz* file.

#### Preprocessing
- *Snakefile* - which contains rules for:
  1. Quality filtering using **AdapterRemoval**
  2. Duplicate read removal using **seqkit rmdup**
  3. Mapping reads against reference genome(s) using **bwa mem**

- Config file *config.yaml*, in which the user may be interested to customise:
  1. Quality filtering - specific adapter sequences, minimum quality, character separating the mate read number.
  2. Mapping reads against reference genome(s) - reference genome(s) path(s), stringent level for mapping and other parameters. 


#### Metagenomics (Individual Assembly & Coassembly)
- *Snakefile* - which contains rules for:
  1. Metagenomic assembly using **metaSpades** or **megahit**
  2. Read mapping to assembly using **bwa mem** 
  3. Contig binning using **Metabat**, **MaxBin** (and **Concoct** for coassembly)
  4. Binner result integration using **DasTool** 
  
- Config file *config.yaml*, in which the user may be interested to customise:
  1. Metagenomic assembly - choose between the mentioned options by writing *megahit* or *spades*
  2. Minimum contig length - minimum bp per contig in final assembly file.

#### Metagenomics Dereplication & Annotation
- *Snakefile* - which contains rules for:
  1. Bin Dereplication using **dRep**
  2. Bin Annotation with **prokka**
  3. Taxonomic classification and phylogenetic inference of bins **GTDB-Tk**

#### Metagenomics Final Statistics
- *Snakefile* - which contains rules for:
  1. Read mapping to MAGs
  2. MAG and contig coverage in reads 



## Exectute Holoflow *.py* workflow launchers
These should be **executed as jobs**, therefore a *.sh* script should be generated which will contain the job itself:

- *.sh* example script for *preprocessing.py* called ***first_job_preprocessing.sh***:
```bash
python full/path/holoflow/preprocessing.py -f full/path/input.txt -d full/path/workdir -c full/path/config.yaml -l full/path/log_file.log -t 40
```

- *job execution* in Computerome2 example:
```bash
 qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e full/path/job_error_file.err -o full/path/job_out_file.out -l nodes=1:ppn=40,mem=180gb,walltime=5:00:00:00 -N JOB_ID full/path/first_job_preprocessing.sh

```
  Note that the job parameters: *ppn*, *nodes*, *memory*, *wall time* ... can and ought to be customised optimally for every job type.





