# holoflow
Bioinformatics pipeline for hologenomics data generation and analysis

Snakemake is a workflow management system which requires from a *Snakefile* and a *config* file. This is a Bioinformatics pipeline for hologenomics data generation and analysis implemented with Snakemake.

## Files and directories
### Main directory

The main *holoflow* directory contains a given number of Python scripts which work as launchers for the different **workflow programs** in the pipeline:

  - *preparegenomes.py*   - Merge all potential reference genomes to sample into a single *.fna* file to be used in preprocessing.py.  
  - *preprocessing.py*    - Data preprocessing from quality to duplicate sequences for further downstream analysis.
  - *metagenomics_IA.py*  - Individual assembly-based assembly and metagenomics binning. 
  
  
These are designed to be called from the command line and require  the following arguments:  

  1. **-f** Input.txt file to *.py* files, which will be used to retrieve fundamental information for the pipeline run.   
  2. **-d** Directory where the pipeline temporary files and directories will be.
  3. **-l** Desired pipeline *log file* path.
  4. **-c** *config* file full path.
  5. **-t** Maximum number of threads to be used by Snakemake.  
  
  
  
#### Input files description
Find *input.txt* file description for every workflow.  
In all cases, columns must be delimited by a simple space and no blank lines should be found in the end of the file.  
Those lines starting by # won't be considered.  
  
##### *preparegenomes.py*

#Genome_ID(nospaces,no-anything) PathGenome NameOutputDB

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
#### Preprocessing
- *Snakefile* - which contains rules for:
  1. Quality filtering using **AdapterRemoval**
  2. Duplicate read removal using **seqkit rmdup**
  3. Mapping reads against reference genome(s) using **bwa mem**

- Config file *config.yaml*, in which the user may be interested to customise:
  1. Quality filtering - specific adapter sequences, minimum quality
  2. Mapping reads against reference genome(s) - reference genome for host and human paths


#### Metagenomics
- *Snakefile* - which contains rules for:
  1. Metagenomic assembly using **metaSpades** or **megahit**
  2. Read mapping to assembly using **bwa mem** ##### UNDER CONSTRUCTION
  3. Contig binning using **Metabat**, **MaxBin** and **Concoct** ##### UNDER CONSTRUCTION
  4. Binner result integration using **DasTool** ##### UNDER CONSTRUCTION
  5. Complementess improvement ##### UNDER CONSTRUCTION
  5. Taxonomic refinement using CAT ##### UNDER CONSTRUCTION
  6. Redundancy refinement ##### UNDER CONSTRUCTION
  7. Dereplication using dRep ##### UNDER CONSTRUCTION
  7. Bin assembly improvement (contig elongation and scaffolding) using SSPACE. ##### UNDER CONSTRUCTION

- Config file *config.yaml*, in which the user may be interested to customise:
  1. Metagenomic assembly - choose between the mentioned options by writing *megahit* or *spades*
  2. Minimum contig length - minimum bp per contig in final assembly file.


## Exectute *holoflow.py*
**The python script should be launched from its containing directory:**
```
python holoflow.py -f ${input} -d ${workdir} -w metagenomics -c ${configfile} -t 40
```
*input*, *workdir* and *configfile* are shell variables which where previously defined in the command line, but the corresponding path to the file can also be directly specified in the python command.


  


  - ***preparegenomes.py*** -


  - ***preprocessing.py*** -


  - ***metagenomics_IA.py*** -




