import argparse
import subprocess
import os
import glob
import sys
import time


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
parser.add_argument('-W', help="rewrite everything", dest="REWRITE", action='store_true')
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
cores=args.threads


    # retrieve current directory
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)

current_time = time.strftime("%m.%d.%y_%H:%M", time.localtime())
if not (args.config_file):
    cpconfigCmd= 'cp '+curr_dir+'/workflows/metagenomics/dereplication/config.yaml '+path+'/'+current_time+'_config.yaml'
    subprocess.Popen(cpconfigCmd,shell=True).wait()

    config = path+'/'+current_time+'_config.yaml'
else:
    config=args.config_file

# If the user does not specify a log file, provide default path
if not (args.log):
    log = os.path.join(path,"Holoflow_dereplication_metagenomics.log")
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
    in_dir = os.path.join(path,"MDR_00-InputBins")

    if not os.path.exists(in_dir): # either because of rewrite or because first time
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Paste desired output file names from input.txt
        group = ''
        output_files=''


        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))

        last_line = lines[-1]

    for line in lines:

        if not (line.startswith('#')):
            dir = line.strip('\n').split(' ') # Create a list of each line

            # the input will be a directory, where all bins for all samples will be contained
            # If Bins from different samples are in different directories, create input Dir
            # and move them all there

            current_input_dir=os.path.dirname(dir[1])
            current_in_files = ''.join(glob.glob(dir[1]+'/*')[1])

            desired_input=(str(in_dir)+'/'+str(dir[0])) # desired input dir path
            if os.path.exists(desired_input):
                desired_in_files = os.listdir(desired_input)

            if args.REWRITE:
                if os.path.basename(current_in_files) in desired_in_files: # the directory has not been yet removed: this group's files already exist in dir
                    rmCmd='rm -rf '+desired_input+''
                    subprocess.Popen(rmCmd,shell=True).wait()
                else:                              # the directory has been  removed already by a previous line in the input file
                    pass

            #if bins not in desired input dir, copy them there
            if not desired_input == current_input_dir:

                if (len(os.listdir(desired_input)) == 0): # if dir exists but empty
                    copyfilesCmd='find  '+dir[1]+' -maxdepth 1 -type f | xargs -I {} ln -s {} '+desired_input+''
                    subprocess.check_call(copyfilesCmd, shell=True)

                if not (os.path.exists(str(desired_input))):
                    copyfilesCmd='mkdir '+desired_input+' && find  '+dir[1]+' -maxdepth 1 -type f | xargs -I {} ln -s {} '+desired_input+''
                    subprocess.check_call(copyfilesCmd, shell=True)

                # write output files

            if not (group == dir[0]): # when the group changes, define output files for previous group#same as last output in Snakefile
                group=str(dir[0])
                final_temp_dir="MDR_03-BinPhylogeny"
                output_files+=(path+"/"+final_temp_dir+"/"+group+"_BAC_Holoflow.gtdbtk_sub.tree ")
                output_files+=(path+"/"+final_temp_dir+"/"+group+"_AR_Holoflow.gtdbtk_sub.tree ")

            if (line == last_line):
                #same as last output in Snakefile
                group=str(dir[0])
                final_temp_dir="MDR_03-BinPhylogeny"
                output_files+=(path+"/"+final_temp_dir+"/"+group+"_BAC_Holoflow.gtdbtk_sub.tree ")
                output_files+=(path+"/"+final_temp_dir+"/"+group+"_AR_Holoflow.gtdbtk_sub.tree ")



    return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/dereplication/Snakefile')

    # Run snakemake
    log_file = open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Metagenomics - Dereplication starting")
    log_file.close()

    mtg_snk_Cmd = 'snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.Popen(mtg_snk_Cmd, shell=True).wait()

    log_file = open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Metagenomics - Dereplication has finished :)")
    log_file.close()

    # Keep temp dirs / remove all
    if args.keep: # If -k, True: keep
        pass
    else: # If not -k, keep only last dir
        exist=list()
        for file in out_files.split(" "):
            exist.append(os.path.isfile(file))

        if all(exist): # all output files exist
            rmCmd='cd '+path+' | grep -v '+final_temp_dir+' | xargs rm -rf && mv '+final_temp_dir+' MDR_Holoflow'
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
