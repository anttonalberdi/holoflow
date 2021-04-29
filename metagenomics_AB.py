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
parser.add_argument('-k', help="keep tmp directories", dest="keep", action='store_true')
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
    config = os.path.join(os.path.abspath(curr_dir),"workflows/metagenomics/individual_binning/config.yaml")
else:
    config=args.config_file
# If the user does not specify a log file, provide default path
if not (args.log):
    log = os.path.join(path,"Holoflow_individualA_metagenomics.log")
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
    in_dir_0 = os.path.join(path,"PPR_03-MappedToReference")

    if not os.path.exists(in_dir_0):
        os.makedirs(in_dir_0)

    with open(in_f,'r') as in_file:
        # Define variables
        output_files=''
        final_temp_dir="MIB_04-BinMerging"
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
                    sample_name=line[0]
                    in_for=line[1]# input for (read1) file
                    in_rev=line[2] # input reverse (read2) file

                    # Define input file
                    in1=in_dir+'/'+sample_name+'_1.fastq'
                    # Check if input files already in desired dir
                    if os.path.isfile(in1) or os.path.isfile(in1+'.gz'):
                        pass
                    else:
                        #If the file is not in the working directory, create soft link in it
                        if os.path.isfile(in_for):
                            if in_for.endswith('.gz'):# if compressed, decompress in standard dir with std ID
                                read1Cmd = 'ln -s '+in_for+' '+in1+'.gz && gunzip -c '+in1+'.gz > '+in1+''
                                subprocess.Popen(read1Cmd, shell=True).wait()
                            else:
                                read1Cmd = 'ln -s '+in_for+' '+in1+''
                                subprocess.Popen(read1Cmd, shell=True).wait()



                    # Define input file
                    in2=in_dir+'/'+sample_name+'_2.fastq'
                    # Check if input files already in desired dir
                    if os.path.isfile(in2) or os.path.isfile(in2+'.gz'):
                        pass
                    else:
                        #If the file is not in the working directory, transfer it
                        if os.path.isfile(in_rev):
                            if in_for.endswith('.gz'):
                                read2Cmd = 'ln -s '+in_rev+' '+in2+'.gz && gunzip -c '+in2+'.gz > '+in2+''
                                subprocess.Popen(read2Cmd, shell=True).wait()
                            else:
                                read2Cmd = 'ln -s '+in_rev+' '+in2+''
                                subprocess.Popen(read2Cmd, shell=True).wait()


                output_files+=(path+"/"+final_temp_dir+"/"+sample_name+"_DASTool_files ")


        else: # the input directory already exists and is full, don't want to create it again, just re-run from last step
            for line in lines:
                ### Skip line if starts with # (comment line)
                if not (line.startswith('#')):

                    line = line.strip('\n').split(' ') # Create a list of each line
                    sample_name=line[0]
                    in_for=line[1]
                    in_rev=line[2]

                output_files+=(path+"/"+final_temp_dir+"/"+sample_name+"_DASTool_files ")



    return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/individual_binning/Snakefile')

    # Run snakemake
    log_file = open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Metagenomics-IndividualBinning starting")
    log_file.close()

    mtg_snk_Cmd = 'snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    log_file = open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Metagenomics-IndividualBinning has finished :)")
    log_file.close()

    # Keep temp dirs / remove all
    if args.keep: # If -k, True: keep
        pass
    else: # If not -k, keep only last dir
        exist=list()
        for file in out_files.split(" "):
            exist.append(os.path.isfile(file))

        if all(exist): # all output files exist
            rmCmd='cd '+path+' | grep -v '+final_temp_dir+' | xargs rm -rf && mv '+final_temp_dir+' MIB_Holoflow'
            subprocess.Popen(rmCmd,shell=True).wait()

        else:   # all expected output files don't exist: keep tmp dirs
            log_file = open(str(log),'a+')
            log_file.write("Looks like something went wrong...\n\t\t The temporal directories have been kept, you should have a look...")
            log_file.close()



###########################
#### Workflows running
###########################
# 2    # Metagenomics workflow
run_metagenomics(in_f, path, config, cores)
