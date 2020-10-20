import argparse
import subprocess
import os
import sys
import ruamel.yaml

###########################
#Argument parsing
###########################
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="input.txt file", dest="input_txt", required=True)
parser.add_argument('-d', help="temp files directory path", dest="work_dir", required=True)
parser.add_argument('-c', help="config file", dest="config_file", required=False)
parser.add_argument('-l', help="pipeline log file", dest="log", required=False)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
cores=args.threads

    # retrieve current directory
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)


if not (args.config_file):
    config = os.path.join(os.path.abspath(curr_dir),"workflows/preprocessing/config.yaml")
else:
    config=args.config_file

if not (args.log):
    log = os.path.join(path,"Holoflow_preprocessing.log")
else:
    log=args.log


    #Append current directory to .yaml config for standalone calling
yaml = ruamel.yaml.YAML()
yaml.explicit_start = True
with open(str(config), 'r') as config_file:
    data = yaml.load(config_file)
    if data == None:
        data = {}

with open(str(config), 'w') as config_file:
    data['holopath'] = str(curr_dir)
    data['logpath'] = str(log)
    dump = yaml.dump(data, config_file)




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
        final_temp_dir="PPR_03-MappedToReference"

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

                    # Add stats and bam output files only once per sample
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+".stats ")
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+"_ref.bam ")

        return output_files



def run_preprocessing(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_preprocessing(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/preprocessing/Snakefile')

    # Run snakemake
    prep_snk_Cmd = 'module unload gcc && module load tools anaconda3/4.4.0 && snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(prep_snk_Cmd, shell=True)
    print("Have a nice run!\n\t\tHOLOFOW Preprocessing starting")


###########################
#### Workflows running
###########################


# 1    # Preprocessing workflow
run_preprocessing(in_f, path, config, cores)
