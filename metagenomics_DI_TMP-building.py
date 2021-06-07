# 02.06.21
import argparse
import subprocess
import os
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
    cpconfigCmd= 'cp '+curr_dir+'/workflows/metagenomics/dietary_analysis/config.yaml '+path+'/'+current_time+'_config.yaml'
    subprocess.Popen(cpconfigCmd,shell=True).wait()

    config = path+'/'+current_time+'_config.yaml'
else:
    config=args.config_file
# If the user does not specify a log file, provide default path
if not (args.log):
    log = os.path.join(path,"Holoflow_DietaryAnalysis_metagenomics.log")
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
    # Find the databases installed by Bent Petersen for annotation of predicted ORFs
    data['db_dir'] = str('/home/projects/ku-cbd/people/nurher/diet_analysis/Diet_DBs')
    data['threads'] = str(cores)
    data['holopath'] = str(curr_dir)
    data['logpath'] = str(log)
    dump = yaml.dump(data, config_file)

###########################
## Functions
###########################



    ###########################
    ###### METAGENOMIC FUNCTIONS

def in_out_dietary_analysis(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    # Define input directory and create it if not exists "00-InputData"
    in_dir = os.path.join(path,"MDI_00-InputData")

    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))

        # Define variables
        output_files=''
        final_temp_dir="MDI_03-Quantify"

    for line in lines:
        ### Skip line if starts with # (comment line)
        if not (line.startswith('#')):

            line = line.strip('\n').split(' ') # Create a list of each line
            group_name=line[0]
            assembly_path=line[1]
            nonmapp_fastq_dir=line[2]

            in_group = in_dir+'/'+group_name
            if os.path.exists(in_group):
                if args.REWRITE:    # if rewrite, remove directory - start from 0
                    rmCmd='rm -rf '+in_group+''
                    subprocess.Popen(rmCmd,shell=True).wait()
                else:               # don't want to rewrite, continue from last rule completed
                    pass

            if not os.path.exists(in_group): # if dir not exists either because of REWRITE or bc first time, DO EVERYTHING
                os.makedirs(in_group)


            # Define output files based on input.txt
            output_files+=path+'/'+final_temp_dir+'/'+group_name+' '

            # Soft link from assembly file
            a_file = in_group+'/'+'group_name.fna'
            if not os.path.isfile(a_file):
                linkAssemblyCmd = 'ln -s '+assembly_path+' '+in_group+'/'+group_name+'.fa'
                subprocess.Popen(linkAssemblyCmd,shell=True).wait()

            # Link .fastq files of non-MAG mapped reads to subdir
            input_nonmapp_dir = in_group+'/'+'mag_unmapped_fastq'

            # Check if input files already in desired dir  -> link fastq of non mapped to MAG reads
            if os.path.exists(input_nonmapp_dir):
                try:    # try to create the link - if the link already exists ... -> TRY/Except is to avoid exception errors
                    mvreadsCmd = 'ln -s '+nonmapp_fastq_dir+'/*notMAGmap*fastq* '+input_nonmapp_dir+''
                    subprocess.Popen(mvreadsCmd, shell=True).wait()
                except: # ... it won't be created, but pass
                    pass
            else:
                mvreadsCmd = 'mkdir '+input_nonmapp_dir+' && ln -s '+nonmapp_fastq_dir+'/*notMAGmap*fastq* '+input_nonmapp_dir+''
                subprocess.Popen(mvreadsCmd, shell=True).wait()

    return output_files


def run_dietary_analysis(in_f, path, config, cores):
    """Run snakemake on shell, wait for it to finish.
    Given flag, decide whether keep only last directory."""

    # Define output names
    out_files = in_out_dietary_analysis(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/dietary_analysis/Snakefile')

    # Run snakemake
    log_file = open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Dietary Analysis starting")
    log_file.close()

    dietary_analysis_snk_Cmd = 'module load tools anaconda3/4.4.0 && snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.Popen(dietary_analysis_snk_Cmd, shell=True).wait()

    log_file = open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Dietary Analysis has finished :)")
    log_file.close()

    # Keep temp dirs / remove all
    if args.keep: # If -k, True: keep
        pass
    else: # If not -k, keep only last dir
        exist=list()
        for file in out_files.split(" "):
            exist.append(os.path.isfile(file))

        if all(exist): # all output files exist
            rmCmd='cd '+path+' | grep -v '+final_temp_dir+' | xargs rm -rf && mv '+final_temp_dir+' MDI_Holoflow'
            subprocess.Popen(rmCmd,shell=True).wait()

        else:   # all expected output files don't exist: keep tmp dirs
            log_file = open(str(log),'a+')
            log_file.write("Looks like something went wrong...\n\t\t The temporal directories have been kept, you should have a look...")
            log_file.close()




###########################
#### Workflows running
###########################


# 1    # Final Stats workflow
run_dietary_analysis(in_f, path, config, cores)
