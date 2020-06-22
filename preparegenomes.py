import argparse
import subprocess
import os
import sys
import re
import ruamel.yaml

###########################
#Argument parsing
###########################
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="input.txt file", dest="input_txt", required=True)
parser.add_argument('-d', help="temp files directory path", dest="work_dir", required=True)
parser.add_argument('-c', help="config file", dest="config_file", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
config=args.config_file
cores=args.threads



    # retrieve current directory
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)

    #Append current directory to .yaml config for standalone calling
yaml = ruamel.yaml.YAML()
yaml.explicit_start = True
with open(str(config), 'r') as config_file:
  data = yaml.load(config_file)

with open(str(config), 'w') as config_file:
  data['holopath'] = str(curr_dir)
  dump = yaml.dump(data, config_file)



###########################
## Functions
###########################

    ###########################
    ###### PREPAREGENOMES FUNCTIONS

def in_out_preparegenomes(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    db_dir = os.path.join(path,"PRG")
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    with open(in_f,'r') as in_file:
        # Paste desired output file names from input.txt
        ref_genomes_IDs=list()
        ref_genomes_paths=list()
        db_ID=''


        lines = in_file.readlines() # Read input.txt lines
        last_file = lines[-1]
        for file in lines:

            if not (file.startswith('#')):
                file = file.strip('\n').split(' ') # Create a list of each line

                # Save IDs for reformat and paths for merging
                ref_genomes_IDs.append(file[0])
                ref_genomes_paths.append(file[1])

                # If all previous genomes to same db, only save db name once
                    # do the merging of the genomes into db
                if (not (re.match(file[2], db_ID))):
                    db_ID = file[2]
                    # call merging function
                    merge_genomes(db_dir,refg_IDs,refg_Paths,db_ID)

                # If ending of lines, and no new db name, also
                    # do the merging of the genomes into db
                if (file == last_file):
                    # call merging function
                    merge_genomes(db_dir,refg_IDs,refg_Paths,db_ID)


> append db path to config ###############################


def merge_genomes(db_dir,refg_IDs,refg_Paths,db_ID):

    for (i in range(len(refg_Paths))):

        genome = refg_Paths[i]
        ID = refg_IDs[i]

        if genome.endswith('.gz'): # uncompress genome for editing
            uncompressCmd='gunzip -c '+genome+' > db_dir###############################'
            subprocess.check_call(uncompressCmd, shell=True)
        else:
            pass

        # edit > genome identifiers
            # find all lines starting with > and add ID_ before all already there
            > save as temp file in db_dir

    # reformat every genome - grep and sed

    # take work dir path and ID and subprocess merge
        > merge all temp files > db_dir/ID.fna
        > remove uncompressed+modified genomes in dir

    return(db_path)  ###############################










def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/individual_assembly/Snakefile')

    # Run snakemake
    mtg_snk_Cmd = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    print("Have a nice run!\n\t\tHOLOFOW Metagenomics starting")


###########################
#### Snakemake pipeline run - load required modules
###########################
load_modulesCmd='module unload gcc/5.1.0 && module load anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)



###########################
#### Workflows running
###########################
# 2    # Metagenomics workflow
run_metagenomics(in_f, path, config, cores)