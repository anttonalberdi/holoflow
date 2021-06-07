import argparse
import subprocess
import glob
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

# If the user does not specify a config file, provide default file in GitHub
current_time = time.strftime("%m.%d.%y_%H:%M", time.localtime())
if not (args.config_file):
    cpconfigCmd= 'cp '+curr_dir+'/workflows/metagenomics/final_stats/config.yaml '+path+'/'+current_time+'_config.yaml'
    subprocess.Popen(cpconfigCmd,shell=True).wait()

    config = path+'/'+current_time+'_config.yaml'
else:
    config=args.config_file

# If the user does not specify a log file, provide default path
if not (args.log):
    log = os.path.join(path,"Holoflow_final_stats.log")
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
    data['KO_DB'] = str('/home/databases/ku-cbd/aalberdi/prokka2kegg/idmapping_KO.tab.gz')
    data['KO_list'] = str(curr_dir+'/workflows/metagenomics/final_stats/KO_list.txt')
    dump = yaml.dump(data, config_file)




###########################
## Functions
###########################



    ###########################
    ###### METAGENOMIC FUNCTIONS

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
        final_temp_dir="MFS_04-KOAbundances"

    for line in lines:
        ### Skip line if starts with # (comment line)
        if not (line.startswith('#')):

            line = line.strip('\n').split(' ') # Create a list of each line
            sample_name=line[0]
            mtg_reads_dir=line[1]
            mtg_files = ''.join(glob.glob(mtg_reads_dir+'/*')[1]) # keep only second metagenomic file
            drep_bins_dir=line[2]
            annot_dir=line[3]

            in_sample = in_dir+'/'+sample_name
            if os.path.exists(in_sample):
                in_mtg_files = os.listdir(in_sample+'/metagenomic_reads') # if the dir already exists, save names of files inside

            if args.REWRITE:    # if rewrite, remove directory
                if os.path.basename(mtg_files) in in_mtg_files: # the directory has not been yet removed: this group's files already exist in dir
                    rmCmd='rm -rf '+in_sample+''
                    subprocess.Popen(rmCmd,shell=True).wait()
                else:                              # the directory has been  removed already by a previous line in the input file
                    pass                           # belonging to the same group, this is the fill-up round

            if not os.path.exists(in_sample): # if dir not exists either because of REWRITE or bc first time, DO EVERYTHING
                os.makedirs(in_sample)
            else:
                pass

            # Define output files based on input.txt
            output_files+=path+'/'+final_temp_dir+'/'+sample_name+' '

            # Define input dir
            in1=in_sample+'/metagenomic_reads'
            # Check if input files already in desired dir
            if os.path.exists(in1):
                try:    # try to create the link - if the link already exists ... -> TRY/Except is to avoid exception errors
                    mvreadsCmd = 'ln -s '+mtg_reads_dir+'/*.fastq* '+in1+''
                    subprocess.Popen(mvreadsCmd, shell=True).wait()
                except: # ... it won't be created, but pass
                    pass
            else:
                mvreadsCmd = 'mkdir '+in1+' && ln -s '+mtg_reads_dir+'/*.fastq* '+in1+''
                subprocess.Popen(mvreadsCmd, shell=True).wait()

# same for the two other directories that have to be created for input

            # Define input dir
            in2=in_sample+'/dereplicated_bins'
            # Check if input files already in desired dir
            if os.path.exists(in2):
                try:
                    mvbinsCmd = 'ln -s '+drep_bins_dir+'/*.fa '+in2+' && cp '+drep_bins_dir+'/../final_bins_Info.csv '+in2+' && cp '+drep_bins_dir+'/../data_tables/Widb.csv '+in2+''
                    subprocess.Popen(mvbinsCmd, shell=True).wait()
                except:
                    pass
            else:
                mvbinsCmd = 'mkdir '+in2+' && ln -s '+drep_bins_dir+'/*.fa '+in2+' && cp '+drep_bins_dir+'/../final_bins_Info.csv '+in2+' && cp '+drep_bins_dir+'/../data_tables/Widb.csv '+in2+''
                subprocess.Popen(mvbinsCmd, shell=True).wait()

            # Define input dir
            in3=in_sample+'/annotation'
            # Check if input files already in desired dir
            if os.path.exists(in3):
                try:
                    mvgffCmd = 'ln -s '+annot_dir+'/*.gff '+in3+''
                    subprocess.Popen(mvgffCmd, shell=True).wait()
                except:
                    pass
            else:
                mvgffCmd = 'mkdir '+in3+' && ln -s '+annot_dir+'/*.gff '+in3+''
                subprocess.Popen(mvgffCmd, shell=True).wait()


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
