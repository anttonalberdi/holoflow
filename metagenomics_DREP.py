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
    config = os.path.join(os.path.abspath(curr_dir),"workflows/metagenomics/dereplication/config.yaml")
else:
    config=args.config_file

if not (args.log):
    log = os.path.join(path,"Holoflow_dereplication_metagenomics.log")
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

    if data['SSPACE']:
        scaffold=True
    else:
        scaffold=False


###########################
## Functions
###########################

    ###########################
    ###### METAGENOMICS FUNCTIONS

def in_out_metagenomics(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    in_dir = os.path.join(path,"MIB_04-BinMerging")
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Paste desired output file names from input.txt
        group = 'empty'
        output_files=''

        if scaffold:
            final_temp_dir="MDRP_03-MAGPhylogenetics"
        if not scaffold:
            final_temp_dir="MDRP_02-MAGPhylogenetics"

        lines = in_file.readlines() # Read input.txt lines
        last_line = lines[-1]
        for line in lines:

            if not (line.startswith('#')):
                dir = line.strip('\n').split(' ') # Create a list of each line

                # the input will be a directory, where all bins for all samples will be contained
                # If Bins from different samples are in different directories, create input Dir
                # and move them all there

                desired_input=(str(in_dir)+'/'+str(dir[0])) # desired input dir path
                current_input_dir=os.path.dirname(dir[1])

                #if bins not in desired input dir, copy them there
                if not desired_input == current_input_dir:
                    if not (os.path.exists(str(desired_input))):
                        os.mkdir(desired_input)
                    else:
                        copyfilesCmd='cp '+dir[1]+' '+desired_input+''
                        subprocess.check_call(copyfilesCmd, shell=True)
                else:
                    pass

                    # write output files
                if group == 'empty': # will only happen on the first round - first group
                    group=str(dir[0])

                elif ((not (group == file[1])) or (line == last_line)): # when the group changes, define output files for previous group
                    #same as last output in Snakefile
                    ####output_files+=?????????(path+"/"+final_temp_dir+"/"+group+" ")
                    group=dir[0] # define new group in case first condition
                    pass


        return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/dereplication/Snakefile')

    # Run snakemake
    mtg_snk_Cmd = 'module unload gcc/5.1.0 && module load tools anaconda3/4.4.0 && snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    print("Have a nice run!\n\t\tHOLOFOW Metagenomics-Dereplication starting")



###########################
#### Workflows running
###########################
# 2    # Metagenomics workflow
run_metagenomics(in_f, path, config, cores)
