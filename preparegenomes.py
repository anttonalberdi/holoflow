import argparse
import subprocess
import os
import glob
import sys
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
    if data == None:
        data = {}

with open(str(config), 'w') as config_file:
    data['holopath'] = str(curr_dir)
    dump = yaml.dump(data, config_file)


###########################
## Functions
###########################

    ###########################
    ###### PREPAREGENOMES FUNCTIONS

def set_up_preparegenomes(path,in_f):
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
        db_paths=''


        lines = in_file.readlines() # Read input.txt lines
        last_file = lines[-1]
        for file in lines:

            if not (file.startswith('#')):
                refg = file.strip('\n').split(' ') # Create a list of each line
                # Save IDs for reformat and paths for merging
                ref_genomes_IDs.append(refg[0])
                ref_genomes_paths.append(refg[1])
                db_ID = refg[2]

                # If all previous genomes to same db, only save db name once
                    # do the merging of the genomes into db
                if not (refg[2] == db_ID):
                    # call merging function
                    db_paths+=''+merge_genomes(db_dir,ref_genomes_IDs,ref_genomes_paths,db_ID)+' '
                    db_ID = refg[2]

                # If ending of lines, and no new db name, also
                    # do the merging of the genomes into db
                if (file == last_file):
                    db_ID = refg[2]
                    # call merging function
                    db_paths+=''+merge_genomes(db_dir,ref_genomes_IDs,ref_genomes_paths,db_ID)+' '

                else:
                    pass

        return(db_paths)






def merge_genomes(db_dir,refg_IDs,refg_Paths,db_ID):

    for i in range(len(refg_Paths)):

        genome = refg_Paths[i]
        ID = refg_IDs[i]

        if genome.endswith('.gz'): # uncompress genome for editing
                                # and save it in db_dir
            uncompressCmd='gunzip -c '+genome+' > '+db_dir+'/'+ID+'.fna'
            subprocess.check_call(uncompressCmd, shell=True)

            # edit ">" genome identifiers
            # find all lines starting with > and add ID_ before all info
            editgenomeCmd='sed -i "s/>/>'+str(ID)+'_/g" '+db_dir+'/'+ID+'.fna'
            subprocess.check_call(editgenomeCmd, shell=True)


        else:
            # move to project dir and edit ">" genome identifiers
            mvgenomeCmd='mv '+genome+' '+db_dir+'/'+ID+'.fna'
            subprocess.check_call(mvgenomeCmd, shell=True)
            editgenomeCmd='sed -i "s/>/>'+str(ID)+'_/g" '+db_dir+'/'+ID+'.fna'
            subprocess.check_call(editgenomeCmd, shell=True)


    # define full db path and merge all reference genomes in it
    db_path = ''+db_dir+'/'+db_ID+'.fna'

        # obtain full paths of all edited genomes to merge
    if os.path.exists(db_path):
        rmCmd='rm '+db_path+''
        subprocess.check_call(rmCmd, shell=True)

    mergeCmd='cd '+db_dir+' && cat *.fna > '+db_path+''
    subprocess.check_call(mergeCmd, shell=True)

    return(db_path)




def run_preparegenomes(in_f, path, config, cores):
    """Run snakemake on shell"""


        # retrieve db_path
    db_paths = set_up_preparegenomes(path,in_f)

        # Append db_path to config for Snakefile running
    yaml = ruamel.yaml.YAML()
    yaml.explicit_start = True
    with open(str(config), 'r') as config_file:
        data = yaml.load(config_file)

    with open(str(config), 'w') as config_file:
        data['DB_path'] = str(db_paths)
        dump = yaml.dump(data, config_file)


    # get output files and Snakefile directory
    out_files = ''+path+'/PRG/ok.txt'
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/preparegenomes/Snakefile')

    # Run snakemake
    prg_snk_Cmd = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(prg_snk_Cmd, shell=True)

    print("Have a nice run!\n\t\tHOLOFOW Prepare genomes starting")





###########################
#### Snakemake pipeline run - load required modules
###########################
load_modulesCmd='module unload gcc/5.1.0 && module load anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)



###########################
#### Workflows running
###########################
# 0    # Preparegenomes workflow
run_preparegenomes(in_f, path, config, cores)
