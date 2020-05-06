# holoflow
Bioinformatics pipeline for hologenomics data generation and analysis

Snakemake is a workflow management system which requires from a *Snakefile* and a *config* file. This is a Bioinformatics pipeline for hologenomics data generation and analysis implemented with Snakemake.

## Files and directories
### Main directory
- *holoflow.py* - which contains the script for the pipeline calling.  
This is designed to be called from the command line, and requires the next arguments:  
  1. **-f** Input file - which will contain three columns delimited by a simple space:  
    a. Sample name. 
    b. Assembly group (If not coassembly this field will be ignored - but it is important that is not omitted when writing the input file). 
    c. Original full path/name of input file/s. 
    d. Final output directory name (*Note it must match the output directory name in the workflow's final Snakefile rule*).  
b. Assembly group (If not coassembly this field will be ignored - but it is important that is not omitted when writing the input file).  
  2. **-d** Project path - directory where pipeline outputs will be stored
  3. **-w** Workflow to be run: preprocessing or metagenomics.
  4. **-config** *config* file full path.
  5.  **-cores** Number of cores to be used by Snakemake.

  
### Workflows - specific directories
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
  3. Contig binning by **Metabat and MaxBin** plus binning refinement by **DasTool** ##### UNDER CONSTRUCTION
  
- Config file *config.yaml*, in which the user may be interested to customise:
  1. Metagenomic assembly - choose between the mentioned options by writing *megahit* or *spades*
  2. Minimum contig length - minimum bp per contig in final assembly file.

  
## Exectute *holoflow.py*
**The python script should be launched from its containing directory:**
```
python holoflow.py -f input.txt -d ${workdir} -w metagenomics -config ${configfile} -cores 40
```
*workdir* and *configfile* are shell variables which where previously defined in the command line, but the corresponding path to the file can also be directly specified in the python command. 



