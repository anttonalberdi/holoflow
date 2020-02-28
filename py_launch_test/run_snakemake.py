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


import argparse
import subprocess
import os
import sys

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="W, dest="inputfile", required=True)
parser.add_argument('-d', help="Project directory path", dest="projectpath", required=True)  ##### ADDED IS NECESSARY PROJECT PATH AT LEAST
parser.add_argument('-w', help="Chosen Workflow", dest="workflow", required=True)
parser.add_argument('-c', help="Config file", dest="configfile", required=True)
args = parser.parse_args()

input_file=args.inputfile
projectpath=args.projectpath
workflow=args.workflow
configfile=args.configfile

# A folder 00-InputData is created in .py and the input files specified in input.txt are moved there
  # and their name changed to the one specified in input.txt's first column.
# Paste desired output and give it to snakemake command


# Create "00-RawData/" directory if not exists
input_dir=os.path.join(projectpath,"00-InputData")
if not os.path.exists(input_dir):
    os.makedirs(input_dir)

    with open(str(input_file),'r') as input_file:
        # Paste desired output file names from input.txt
        read = 0
        output_files=''
        for file in input_file:
            file = file.split()
            read+=1
            output_files+=(projectpath+"/"+file[2]+"/"+file[0]+"_"+str(read)+".fastq ")   ####### should be independent from.fastq (TRY UNTIL 04 MAP HUMAN)
            if read == 2:
                read=0

            #Move files to new dir "00-RawData/" and change file names for 1st column in input.txt
            filename=file[1]
            copyfilesCmd='cp '+filename+' '+input_dir+''
            subprocess.check_call(copyfilesCmd, shell=True)

            new_name=file[0]
            renamefilesCmd='cd '+input_dir+' && mv '+filename+' '+new_name+''


# Snakemake pipeline run
load_modulesCmd='module unload gcc/5.1.0 && module load anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)

###########################
######## WORKFLOWS ########
###########################

    # Preprocessing workflow
if workflow == "preprocessing":
    snakemakeCmd = 'xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` '+projectpath+'/snakemake.log -l nodes=1:ppn=28,mem=100gb,walltime=0:06:00:00 -N holoflow_preprocessing -de snakemake -s preprocessing/Snakefile '+output_files+' --configfile '+configfile+''
    subprocess.check_call(snakemakeCmd, shell=True)



    # Metagenomics workflow
if workflow == "metagenomics":

    prep = input("Input files for holoflow/metagenomics are fastq. Is your data preprocessed? write[y/n]")

    if prep == 'n':
        prep2 = input("Would you like to process it before running holoflow/metagenomics with holoflow/preprocessing? write[y/n]")

        if prep2 == 'n':
            print("You should come back when your data is preprocessed. See you soon :)")
        if prep2 == 'y':
            snakemakeCmd = 'xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` '+projectpath+'/snakemake.log -l nodes=1:ppn=28,mem=100gb,walltime=0:06:00:00 -N holoflow_metagenomics -de snakemake -s metagenomics/prep_and_metagenomics/Snakefile '+output_files+' --configfile '+configfile+''
            subprocess.check_call(snakemakeCmd, shell=True)

    if prep == 'y':
        print("Great! Have a nice running!\n\t\tHOLOFOW Metagenomics starting")
        snakemakeCmd = 'xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` '+projectpath+'/snakemake.log -l nodes=1:ppn=28,mem=100gb,walltime=0:06:00:00 -N holoflow_metagenomics -de snakemake -s metagenomics/Snakefile '+output_files+' --configfile '+configfile+''
        subprocess.check_call(snakemakeCmd, shell=True)


    # Genomics workflow
