# holoflow
Bioinformatics pipeline for hologenomics data generation and analysis

Snakemake is a workflow management system which requires from a *Snakefile* and a *config* file. This is a Bioinformatics pipeline for hologenomics data generation and analysis implemented with Snakemake.

## Files and directories
### Main directory
- *run_snakemake.py* - which contains the script for the pipeline calling.  
This is designed to be called from the command line, and requires the next arguments:
  1. **-f** Input file - which will contain three columns delimited by **\t**:
    a. Sample name
    b. Original full path/name of input file
    c. Final output directory name (*Note it must match the output directory name in the desired final Snakefile rule*)
  2. **-d** Project path - where pipeline outputs will be stored
  2. **-w** Workflow to be run
  2. **-c** *config* file full path 

  
### Workflow-specific directories

#### Metagenomics
- *Snakefile* - which contains rules for:
  1. Quality filtering using **AdapterRemoval** 
  2. Duplicate read removal using **seqkit rmdup**
  3. Mapping reads against reference genome(s) using **bwa mem**
  4. Metagenomic assembly using **metaSpades** or **megahit**
  5. Read mapping to assembly using **bwa mem**
  6. Contig binning by **Metabat and MaxBin** plus binning refinement by **DasTool**
  
- Config file *config.yaml*, in which the user may be interested to customise:
  1. Quality filtering - specific adapter sequences, min quality
  2. Mapping reads against reference genome(s) - reference host and human genome paths
  3. Metagenomic assembly - choose between the mentioned options by writing *megahit* or *spades*
  
  ...among others. 
  
## Exectute *run_snakemake.py*
In case the python script is runned from the directory which contains it:
```
python run_snakemake.py -f input.txt -d ${workdir} -w metagenomics -c ${configfile}
```
*workdir* and *configfile* are shell variables which where previously defined in the terminal, but the corresponding path to the file can also be directly specified in the python command. 

############
module unload gcc/5.1.0
module load anaconda3/4.4.0

snakemake -s Snakefile -n -r ${workdir}/02-DuplicatesRemoved/H2A_1.fastq ${workdir}/02-DuplicatesRemoved/H2A_2.fastq
  
