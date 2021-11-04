import argparse
import subprocess
import os
import sys

###########################
#Argument parsing
###########################
# Gather input files and variables from command line
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="input.txt file", dest="input_txt", required=True)
parser.add_argument('-d', help="temp files directory path", dest="work_dir", required=True)
parser.add_argument('-c', help="config file", dest="config_file", required=False)
parser.add_argument('-l', help="pipeline log file", dest="log", required=False)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-N', help="JOB ID", dest="job", required=True)
parser.add_argument('-W', help="rewrite everything", dest="REWRITE", action='store_true')
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
cores=args.threads
job=args.job

    # retrieve current directory
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)

# If the user does not specify a config file, provide default file in GitHub
if not (args.config_file):
    cpconfigCmd= 'cp '+curr_dir+'/workflows/metagenomics/assembly_based/config.yaml '+path+'/'+job+'_config.yaml'
    subprocess.Popen(cpconfigCmd,shell=True).wait()

    config = path+'/'+job+'_config.yaml'
else:
    config=args.config_file
# If the user does not specify a log file, provide default path
if not (args.log):
    log = os.path.join(path,"Holoflow_AssemblyBased_metagenomics.log")
else:
    log=args.log


    # Load dependencies
loaddepCmd='module unload gcc && module load tools anaconda3/4.4.0'
subprocess.Popen(loaddepCmd,shell=True).wait()

    #Append current directory to .yaml config for standalone calling
    # see preprocessing.py for verbose description
import ruamel.yaml
yaml = ruamel.yaml.YAML()
yaml.explicit_start = True
with open(str(config), 'r') as config_file:
    data = yaml.load(config_file)
    if data == None:
        data = {}

with open(str(config), 'w') as config_file:
    # config file to be loaded by DRAM to find the databases installed by Bent Petersen
    data['DRAM_config'] = str('/home/databases/ku-cbd/DRAM/20210705/20210705.dram.config')
    # provided conda environment file by DRAM developers  -> run DRAM in conda env - module way tries to modify internal paths
    data['conda_env_file'] = str(curr_dir+'/workflows/metagenomics/assembly_based/environment.yaml')
    data['threads'] = str(cores)
    data['holopath'] = str(curr_dir)
    data['logpath'] = str(log)
    dump = yaml.dump(data, config_file)


###########################
## Functions
###########################

    ###########################
    ###### METAGENOMICS FUNCTIONS

def in_out_metagenomics(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    in_dir_0 = os.path.join(path,"MAB_00-InputData") # general path

    if not os.path.exists(in_dir_0):
        os.makedirs(in_dir_0)

    with open(in_f,'r') as in_file:
        # Define variables
        output_files=''
        final_temp_dir="MAB_01-Annotation"
        all_lines = in_file.readlines() # Read input.txt lines

        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))


    if os.path.exists(in_dir_0):  # Already run for: same job (wants to continue/Rewrite), for another job
        # Define specific job dir
        in_dir=in_dir_0+'/'+job
        # Define specific job final output dir - for snakemake (needs output files)
        final_temp_dir=final_temp_dir+'/'+job

        # If user wants to remove previous runs' data and run from scratch
        if args.REWRITE:
            if os.path.exists(in_dir):
                rmCmd='rm -rf '+in_dir+''
                subprocess.Popen(rmCmd,shell=True).wait()

        if not os.path.exists(in_dir): # if specific job input directory does not exist
            os.makedirs(in_dir)

        else: # already exists and don't want to rewrite, then pass
            pass

        # If directory is empty, do all - otherwise, just save output names
        if len(os.listdir(in_dir) ) == 0:

            for line in lines:# for line in lines in input file, do:
                ### Skip line if starts with # (comment line)
                if not (line.startswith('#')):

                    line = line.strip('\n').split(' ') # Create a list of each line
                    assembly_id=line[0]
                    assembly_path=line[1]# input for (read1) file

                    # Define input file
                    in1=in_dir+'/'+assembly_id+'.fa'
                    # Check if input files already in desired dir
                    if os.path.isfile(in1) or os.path.isfile(in1+'.gz'):
                        pass
                    else:
                        #If the file is not in the working directory, create soft link in it
                        if os.path.isfile(assembly_path):
                            if assembly_path.endswith('.gz'):# if compressed, decompress in standard dir with std ID
                                read1Cmd = 'ln -s '+assembly_path+' '+in1+'.gz && gunzip -c '+in1+'.gz > '+in1+''
                                subprocess.Popen(read1Cmd, shell=True).wait()
                            else:
                                read1Cmd = 'ln -s '+assembly_path+' '+in1+''
                                subprocess.Popen(read1Cmd, shell=True).wait()


                output_files+=(path+"/"+final_temp_dir+"/"+assembly_id+" ")


        else: # the input directory already exists and is full, don't want to create it again, just re-run from last step
            for line in lines:
                ### Skip line if starts with # (comment line)
                if not (line.startswith('#')):

                    line = line.strip('\n').split(' ') # Create a list of each line
                    assembly_id=line[0]
                    assembly_path=line[1]

                output_files+=(path+"/"+final_temp_dir+"/"+assembly_id+" ")



    return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)

    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/assembly_based/Snakefile')

    # Run snakemake
    log_file = open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Metagenomics-AssemblyBased starting")
    log_file.close()

    mtg_snk_Cmd = 'snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    log_file = open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Metagenomics-AssemblyBased has finished :)")
    log_file.close()


###########################
#### Workflows running
###########################
# 2    # Metagenomics workflow
run_metagenomics(in_f, path, config, cores)
