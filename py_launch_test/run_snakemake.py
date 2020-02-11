# 1. PYTHON SCRIPT to launch SNAKEMAKE
  # PARSE:
    # -f input.txt file path
    # -w workflow (metagenomics, preprocessing...)
      # create conditions, if workflow object == metagenomics, then call the specific
      # Snakemake file and its path in the holoflow folder and run it
    # -c config.yaml file path
  # Create the option for intermediate folders to be deleted when final output is obtained.
  # A folder 00-RAWDATA is created in .py and the input files specified in input.txt are moved there
    # and their name changed to the one specified in input.txt's first column.
  # In snakemake command, input dir is specified by --config inputdir="bla/bla"
  # Paste output and give it to snakemake command

    # Input.txt that contains three columns:
      # - new name of sample (twice_1,_2) 5001
      # - input path (together with next?)
      # - full name in input directory
      # - desired FINAL output dir (the one to be specified in snakemake command!)
        # ------- In python script open this file and paste (outputdir/newname and give it to Snakemake as output files)

# 2. CONFIG.yaml
  #Change in config.yaml inputdir for path_dir !!!!!!!!!!!!!!!!!

# 3. SNAKEFILE
  # Modify first rule's input for 00-RAWDATA/file (check output info)

import argparse
import subprocess
import os

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="Input file", dest="inputfile", required=True)
parser.add_argument('-d', help="Project directory", dest="projectpath", required=True)  ##### ADDED IS NECESSARY PROJECT PATH AT LEAST 
parser.add_argument('-w', help="Chosen Workflow", dest="workflow", required=True)
parser.add_argument('-c', help="Config file", dest="configfile", required=True)

input=args.inputfile
projectpath=args.projectpath
workflow=args.workflow
config_file=args.configfile

# A folder 00-RAWDATA is created in .py and the input files specified in input.txt are moved there
  # and their name changed to the one specified in input.txt's first column.
# Paste output and give it to snakemake command


# Create "00-RawData/" directory if not exists
input_dir=os.join(projectpath,"00-RawData")
if not os.path.exists(input_dir):
    os.makedirs(input_dir)

#Move files to new dir "00-RawData/" and change names of files for 1st column name
for filename in os.listdir('.'):
    copyfilesCmd='cp '+filename+' '+input_dir+''
    subprocess.check_call(copyfilesCmd, shell=True)


# Paste new name to rename files
    renamefilesCmd='mv '+filename+' '+NEW NAME TO BE PASTED YET+' '



# Paste desired output files from input.txt
with open("input.txt",'r') as input_file:
    read = 0
    output_files=[]
    for i in range(len(input_file)):
        read+=1
        output_files.append(input_file[i][2]+"/"+input_file[i][0]+"_"+str(read)+".fastq")   ####### IS THIS CORRECT? SPECIFY .fastq
        if read == 2:
            read=0



# Snakemake pipeline run
load_modulesCmd='module unload gcc/5.1.0 && module load anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)


    # Metagenomics workflow
if workflow == "metagenomics":
    snakemakeCmd = 'xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` ####-e ${workdir}/snakemake.log#### -l nodes=1:ppn=28,mem=100gb,walltime=0:06:00:00 -N holoflow_test -de snakemake -s metagenomics/Snakefile '+output_files+''
    subprocess.check_call(snakemakeCmd, shell=True)

 ## WHAT TO DO WITH -e PATH/blabla.log PATH?
