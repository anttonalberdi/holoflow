# holoflow
Bioinformatics pipeline for hologenomics data generation and analysis

Snakemake is a workflow management system which requires from a *Snakefile* and a *config* file. This is a Bioinformatics pipeline for hologenomics data generation and analysis implemented with Snakemake.

## Files and directories
### Main directory
- *holoflow.py* - which contains the script for the pipeline calling.  
This is designed to be called from the command line, and requires the following arguments:  
  1. **-f** Input.txt file to holoflow.py, which will be used to retrieve fundamental information for the pipeline run. It must contain three columns delimited by a simple space:  
    a. Sample name.  
    b. Assembly group (If not coassembly this field will be ignored - but it is important that is not omitted when writing the input file).  
    c. Original full path/name of input file/s.  
    d. Final output directory name (*Note it must match the output directory name in the workflow's final Snakefile rule*).    
  2. **-d** Directory where the pipeline temporary files and directories will be.
  3. **-w** Workflow to be run: preprocessing or metagenomics.
  4. **-c** *config* file full path.
  5. **-t** Maximum number of threads to be used by Snakemake.

#### Example of input file
|   |   |   |
| --- | --- | --- |
| Sample1 | Group1 | /home/Sample1_1.fq;/home/Sample1_2.fq |
| Sample2 | Group1 | /home/Sample2_1.fq;/home/Sample1_2.fq |
| Sample3 | Group2 | /home/Sample3_1.fq;/home/Sample3_2.fq |
| Samplen | Groupn | /home/Samplen_1.fq;/home/Samplen_2.fq |
  
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



