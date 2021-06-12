import argparse
import subprocess
import os
import glob
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
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
cores=args.threads

    # retrieve current directory
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)

# If the user does not specify a config file, provide default file in GitHub
# If the user does not specify a config file, provide default file in GitHub
if not (args.config_file):
    cpconfigCmd= 'cp '+curr_dir+'/workflows/preparegenomes/config.yaml '+path+'/config.yaml'
    subprocess.Popen(cpconfigCmd,shell=True).wait()

    config = path+'/config.yaml'
else:
    config=args.config_file

# If the user does not specify a log file, provide default path
if not (args.log):
    log = os.path.join(path,"Holoflow_preparegenomes.log")
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
        output_files=''


        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))

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
                    db_paths+=''+merge_genomes(ref_genomes_IDs,ref_genomes_paths,db_ID)+' '
                    output_files+=''+path+'/'+db_ID+'.tar.gz'
                    db_ID = refg[2]
                    ref_genomes_IDs=list()
                    ref_genomes_paths=list()


                # If ending of lines, and no new db name, also
                    # do the merging of the genomes into db
                if (file == last_file):
                    db_ID = refg[2]
                    # call merging function
                    db_paths+=''+merge_genomes(ref_genomes_IDs,ref_genomes_paths,db_ID)+' '
                    output_files+=''+path+'/'+db_ID+'.tar.gz'

                else:
                    pass

        return[db_paths,output_files]






def merge_genomes(refg_IDs,refg_Paths,db_ID):

    db_dir = os.path.join(path,"PRG")

    if not (os.path.exists(str(''+db_dir+'/'+db_ID+'.fna'))):
        for i in range(len(refg_Paths)):

            genome = refg_Paths[i]
            ID = refg_IDs[i]


            if not (os.path.exists(str(''+db_dir+'/'+ID+'.fna'))):
                if genome.endswith('.gz'): # uncompress genome for editing
                                        # and save it in db_dir

                    uncompressCmd='ln -s '+genome+' '+db_dir+'/'+ID+'.fna.gz && gunzip -c '+db_dir+'/'+ID+'.fna.gz > '+db_dir+'/'+ID+'.fna'
                    subprocess.check_call(uncompressCmd, shell=True)

                    # edit ">" genome identifiers
                    # find all lines starting with > and add ID_ before all info
                    editgenomeCmd='sed -i "s/>/>'+str(ID)+'_/g" '+db_dir+'/'+ID+'.fna'
                    subprocess.check_call(editgenomeCmd, shell=True)


                else:
                    # move to project dir and edit ">" genome identifiers
                    mvgenomeCmd='ln -s '+genome+' '+db_dir+'/'+ID+'.fna'
                    subprocess.check_call(mvgenomeCmd, shell=True)
                    editgenomeCmd='sed -i "s/>/>'+str(ID)+'_/g" '+db_dir+'/'+ID+'.fna'
                    subprocess.check_call(editgenomeCmd, shell=True)


        # define full db path and merge all reference genomes in it
        db_path = ''+db_dir+'/'+db_ID+'.fna'

        # obtain full paths of all edited genomes to merge
        mergeCmd='cd '+db_dir+' && cat *.fna > '+db_path+''
        subprocess.check_call(mergeCmd, shell=True)

        # remove all individual genomes if more than one
        if os.path.exists(db_dir+"/"+ID+".fna"):
            rmCmd='cd '+db_dir+' && ls | grep -v "'+db_ID+'*" | xargs rm'
            subprocess.check_call(rmCmd, shell=True)
        else:
            pass

    else: # the db file already exists
        # define full db path and merge all reference genomes in it
        db_path = ''+db_dir+'/'+db_ID+'.fna'


    return(db_path)




def run_preparegenomes(in_f, path, config, cores):
    """Run snakemake on shell"""


        # retrieve db_path
    path_out = set_up_preparegenomes(path,in_f)

        # Append db_path to config for Snakefile running
    yaml = ruamel.yaml.YAML()
    yaml.explicit_start = True
    with open(str(config), 'r') as config_file:
        data = yaml.load(config_file)

    with open(str(config), 'w') as config_file:
        data['DB_path'] = str(path_out[0]).strip()
        dump = yaml.dump(data, config_file)


    # get output files and Snakefile directory
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/preparegenomes/Snakefile')

    # Run snakemake
    log_file = open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Preparegenomes starting")
    log_file.close()

    prg_snk_Cmd = 'snakemake -s '+path_snkf+' -k '+path_out[1]+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(prg_snk_Cmd, shell=True)

    log_file = open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Preparegenomes has finished :)")
    log_file.close()


    #Check how the run went

    for file in path_out.split(" "):
        exist.append(os.path.isfile(file))

    if not all(exist): # all output files exist

        log_file = open(str(log),'a+')
        log_file.write("Looks like something went wrong...\n\t\t")
        log_file.close()





###########################
#### Workflows running
###########################
# 0    # Preparegenomes workflow
run_preparegenomes(in_f, path, config, cores)
