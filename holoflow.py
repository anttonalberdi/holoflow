import argparse
import subprocess
import os
import sys

###########################
#Argument parsing
###########################
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="input.txt file", dest="input_txt", required=True)
parser.add_argument('-d', help="temp files directory path", dest="work_dir", required=True)
parser.add_argument('-w', help="chosen workflow", dest="workflow", required=True)
parser.add_argument('-c', help="config file", dest="config_file", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
workflow=args.workflow
config=args.config_file
cores=args.threads


###########################
## Functions
###########################

    ###########################
    ###### PREPROCESSING FUNCTIONS

def in_out_preprocessing(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    # Define input directory and create it if not exists "00-InputData"
    in_dir = os.path.join(path,"PPR_00-InputData")
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Generate desired output file names from input.txt
        read = 0
        output_files=''
        final_temp_dir="PPR_04-MappedToHuman"

        lines = in_file.readlines() # Read input.txt lines
        for file in lines:

            if not (file.startswith('#')):
                file = file.strip('\n').split(' ') # Create a list of each line

                read+=1     # every sample will have two reads, keep the name of the file but change the read
                # Add an output file based on input.txt info to a list for Snakemake command
                output_files+=(path+"/"+final_temp_dir+"/"+file[0]+"_"+str(read)+".fastq ")

                # Move files to new dir "00-InputData" and change file names for 1st column in input.txt
                #   if the current input file names do not match the designed ones in input.txt
                filename=file[2]       # current input file path and name
                desired_filename='"'+in_dir+'/'+file[0]+'_'+str(read)+'.fastq"'  # desired input file path and name specified in input.txt

                if not ((filename == desired_filename) and (os.path.exists(str(desired_filename)))):
                    if filename.endswith('.gz'):    # uncompress input file if necessary
                        uncompressCmd='gunzip -c '+filename+' > '+desired_filename+''
                        subprocess.check_call(uncompressCmd, shell=True)
                    else:                           # else just move the input file to "00-InputData" with the new name
                        copyfilesCmd='cp '+filename+' '+desired_filename+''
                        subprocess.check_call(copyfilesCmd, shell=True)


                if read == 2:
                    read=0  # two read files for one sample finished, new sample

                    # Add stats output file only once per sample
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+".stats ")

        return output_files



def run_preprocessing(in_f, path, config, cores):
    """Run snakemake on shell"""
    
    # Define output names
    out_files = in_out_preprocessing(path,in_f)
    path_snkf = os.path.join(curr_dir,'workflows/preprocessing/Snakefile')

    # Run snakemake
    prep_snk_Cmd = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(prep_snk_Cmd, shell=True)
    print("Have a nice run!\n\t\tHOLOFOW Preprocessing starting")






    ###########################
    ###### METAGENOMICS FUNCTIONS

def in_out_metagenomics(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    in_dir = os.path.join(path,"PPR_04-MappedToHuman")
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Paste desired output file names from input.txt
        read = 0
        output_files=''
        final_temp_dir="MIA_03-Binning"

        lines = in_file.readlines() # Read input.txt lines
        for file in lines:

            if not (file.startswith('#')):
                file = file.strip('\n').split(' ') # Create a list of each line

                read+=1 # every sample will have two reads, keep the name of the file but change the read

                # Add an output file based on input.txt info to a list for Snakemake command
                output_files+=(path+"/"+final_temp_dir+"/"+file[0]+"_dastool/"+file[0])


                # Move files to new dir "PPR_04-MappedToHuman/" and change file names for 1st column in input.txt
                #   if the current input file names do not match the designed ones in input.txt
                filename=file[2]       # current input file path and name
                desired_filename='"'+in_dir+'/'+file[0]+'_'+str(read)+'.fastq"' # desired input file path and name specified in input.txt

                if not ((filename == desired_filename) and (os.path.exists(str(desired_filename)))):
                    if filename.endswith('.gz'): # uncompress input file if necessary
                        uncompressCmd='gunzip -c '+filename+' > '+desired_filename+''
                        subprocess.check_call(uncompressCmd, shell=True)

                    else:   # else just move the input file to "00-InputData" with the new name
                        copyfilesCmd='cp '+filename+' '+desired_filename+''
                        subprocess.check_call(copyfilesCmd, shell=True)


                if read == 2: # two read files for one sample finished, new sample
                    read=0
                    # Add stats output file only once per sample
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+".stats ")

        return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    path_snkf = os.path.join(curr_dir,'workflows/metagenomics/individual_assembly/Snakefile')

    # Run snakemake
    mtg_snk_Cmd = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    print("Have a nice run!\n\t\tHOLOFOW Metagenomics starting")






    ###########################
    ###### PREPROCESSING AND METAGENOMICS FUNCTIONS

# def run_prepandmet(prepin_f, metin_f, path, prepconfig, metconfig, cores):
#     """Run both preprocessing and metagenomics Snakefiles on shell"""
#
#     # Define output names
#     out_files = in_out_metagenomics(path,in_f)
#


###########################
#### Snakemake pipeline run - load required modules
###########################
load_modulesCmd='module unload gcc/5.1.0 && module load anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)



###########################
#### Workflows
###########################

# 0    # Prepare genomes workflow



# 1    # Preprocessing workflow
if workflow == "preprocessing":
    run_preprocessing(in_f, path, config, cores)


# 2    # Metagenomics workflow

if workflow == "metagenomics": # DATA HAS TO BE PREPROCESSED!
    run_metagenomics(in_f, path, config, cores)


# 3    # Genomics workflow
