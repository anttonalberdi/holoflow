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
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
cores=args.threads

    # retrieve current directory
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)


if not (args.config_file):
    config = os.path.join(os.path.abspath(curr_dir),"workflows/metagenomics/final_stats/config.yaml")
else:
    config=args.config_file

if not (args.log):
    log = os.path.join(path,"Holoflow_final_stats.log")
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
    data['holopath'] = str(curr_dir)
    data['logpath'] = str(log)
    dump = yaml.dump(data, config_file)




###########################
## Functions
###########################



    ###########################
    ###### PREPROCESSING FUNCTIONS

def in_out_final_stats(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    # Define input directory and create it if not exists "00-InputData"
    in_dir = os.path.join(path,"MFS_00-InputData")

    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))

        # Define variables
        output_files=''
        final_temp_dir="MFS_02-MAGCoverage"

        for line in lines:
            ### Skip line if starts with # (comment line)
            if not (line.startswith('#')):

                line = line.strip('\n').split(' ') # Create a list of each line
                sample_name=line[0]
                mtg_reads_dir=line[1]
                drep_bins_dir=line[2]

                # Define output files based on input.txt
                output_files+=path+'/'+final_temp_dir+'/'+sample_name+'/'+sample_name+'.coverage_byMAG.txt '

                # Define input dir
                in1=in_dir+'/'+sample_name+'/metagenomic_reads'
                # Check if input files already in desired dir
                if os.path.exists(in1):
                    pass
                else:
                    mvreadsCmd = 'cd '+mtg_reads_dir+' && cp *.fastq '+in1+''
                    subprocess.Popen(mvreadsCmd, shell=True).wait()


                # Define input dir
                in2=in_dir+'/'+sample_name+'/dereplicated_bins'
                # Check if input files already in desired dir
                if os.path.exists(in2):
                    pass
                else:
                    mvbinsCmd = 'cd '+drep_bins_dir+' && cp *.fa '+in2+''
                    subprocess.Popen(mvbinsCmd, shell=True).wait()


        return output_files



def run_final_stats(in_f, path, config, cores):
    """Run snakemake on shell, wait for it to finish.
    Given flag, decide whether keep only last directory."""

    # Define output names
    out_files = in_out_final_stats(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/final_stats/Snakefile')

    # Run snakemake
    log_file = open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Final Stats starting")
    log_file.close()

    final_stats_snk_Cmd = 'module load tools anaconda3/4.4.0 && snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.Popen(final_stats_snk_Cmd, shell=True).wait()

    log_file = open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Final Stats has finished :)")
    log_file.close()

    # Keep temp dirs / remove all
    if args.keep: # If -k, True: keep
        pass
    else: # If not -k, keep only last dir
        exist=list()
        for file in out_files.split(" "):
            exist.append(os.path.isfile(file))

        if all(exist): # all output files exist
            rmCmd='cd '+path+' | grep -v '+final_temp_dir+' | xargs rm -rf && mv '+final_temp_dir+' MFS_Holoflow'
            subprocess.Popen(rmCmd,shell=True).wait()

        else:   # all expected output files don't exist: keep tmp dirs
            log_file = open(str(log),'a+')
            log_file.write("Looks like something went wrong...\n\t\t The temporal directories have been kept, you should have a look...")
            log_file.close()




###########################
#### Workflows running
###########################


# 1    # Final Stats workflow
run_final_stats(in_f, path, config, cores)