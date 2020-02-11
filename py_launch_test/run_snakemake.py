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
  # In snakemake command, input dir is specified by --config input_dir="bla/bla"
  # Paste output and give it to snakemake command

    # Input.txt that contains three columns:
      # - new name of sample (twice_1,_2) 5001
      # - input path (together with next?)
      # - full name in input directory
      # - desired FINAL output dir (the one to be specified in snakemake command!)
        # ------- In python script open this file and paste (outputdir/newname and give it to Snakemake as output files)

# 2. CONFIG.yaml
  #Change in config.yaml input_dir for path_dir !!!!!!!!!!!!!!!!!

# 3. SNAKEFILE
  # Modify first rule's input for 00-RAWDATA/file (check output info)

import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="Input file", dest="inputfile", required=True)
parser.add_argument('-w', help="Chosen Workflow", dest="workflow", required=True)
parser.add_argument('-c', help="Config file", dest="configfile", required=True)

input=args.inputfile
workflow=args.workflow
config_file=args.configfile
