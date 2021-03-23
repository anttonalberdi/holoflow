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
parser.add_argument('-c', help="config file", dest="config_file", required=False)
parser.add_argument('-k', help="keep tmp directories", dest="keep", action='store_true')
parser.add_argument('-l', help="pipeline log file", dest="log", required=False)
parser.add_argument('-t', help="threads", dest="threads", required=True)
parser.add_argument('-R', help="threads", dest="RERUN", action='store_true')
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


    # Load dependencies
loaddepCmd='module unload gcc && module load tools anaconda3/4.4.0'
subprocess.Popen(loaddepCmd,shell=True).wait()


    #Append current directory to .yaml config for standalone calling
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

    if not os.path.exists(in_dir):
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

        if not args.RERUN: # RE RUN FROM SCRATCH

            if os.path.exists(in_dir):
                rmCmd='rm -rf '+in_dir+''
                subprocess.Popen(rmCmd,shell=True).wait()
                os.makedirs(in_dir)

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
                            copyfilesCmd='mkdir '+desired_input+' && find  '+dir[1]+' -maxdepth 1 -type f | xargs -I {} ln -s {} '+desired_input+''
                            subprocess.check_call(copyfilesCmd, shell=True)

                        if (os.path.exists(str(desired_input))):
                            try:
                                copyfilesCmd='find  '+dir[1]+' -maxdepth 1 -type f | xargs -I {} ln -s {} '+desired_input+''
                                subprocess.check_call(copyfilesCmd, shell=True)
                            except:
                                pass

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



        if args.RERUN: ## RERUN FROM LAST RUN RULE

            for line in lines:
                if not (line.startswith('#')):
                    dir = line.strip('\n').split(' ') # Create a list of each line

                    # the input will be a directory, where all bins for all samples will be contained
                    # If Bins from different samples are in different directories, create input Dir
                    # and move them all there

                    desired_input=(str(in_dir)+'/'+str(dir[0])) # desired input dir path
                    current_input_dir=os.path.dirname(dir[1])

                    if (not (group == dir[0])): # when the group changes, define output files for previous group
                        #same as last output in Snakefile
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
    #subprocess.Popen(mtg_snk_Cmd, shell=True).wait()

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
