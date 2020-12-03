import argparse
import subprocess
import os
import glob
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
    config = os.path.join(os.path.abspath(curr_dir),"workflows/metagenomics/coassembly_binning/config.yaml")
else:
    config=args.config_file

if not (args.log):
    log = os.path.join(path,"Holoflow_coassembly_metagenomics.log")
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
    ###### METAGENOMICS FUNCTIONS

def in_out_metagenomics(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    in_dir = os.path.join(path,"PPR_03-MappedToReference")

    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Define variables
        coa_group = False
        coa1_filename=''
        coa2_filename=''
        read1_files=''
        read2_files=''
        output_files=''
        final_temp_dir="MCB_04-BinMerging"

        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))
        last_line = lines[-1]

        for line in lines:

            if not (line.startswith('#')):
                line = line.strip('\n').split(' ') # Create a list of each line
                sample=str(line[0])      # sample ID

                if (coa_group == line[1]) or not coa_group: # If sample in same coa group or not defined yet
                    read1_files+=line[2]+' '
                    read2_files+=line[3]+' '
                    coa_group=line[1]


                if coa_group and not (coa_group == line[1]): # When the coa group is defined and changes, define output files for previous group and finish input

                    # Finish last coa group
                    coa1_filename=(str(in_dir)+'/'+str(coa_group)+'_1.fastq')
                    coa2_filename=(str(in_dir)+'/'+str(coa_group)+'_2.fastq')

                    # Define Snakemake output files
                    output_files+=(path+"/"+final_temp_dir+"/"+coa_group+"_DASTool_bins ")

                    # Define Snakemake input files
                    if not (os.path.exists(coa1_filename) and os.path.exists(coa2_filename)):
                            # merge all .fastq for coassembly
                        merge1Cmd='cat '+read1_files+' > '+coa1_filename+''
                        subprocess.check_call(merge1Cmd, shell=True)

                        merge2Cmd='cat '+read2_files+' > '+coa2_filename+''
                        subprocess.check_call(merge2Cmd, shell=True)

                    else:
                        pass

                    # Define new coa group
                    coa_group=line[1]
                    read1_files+=line[2]+' '
                    read2_files+=line[3]+' '


                if line == last_line:
                    # Finish last coa group
                    coa1_filename=(str(in_dir)+'/'+str(coa_group)+'_1.fastq')
                    coa2_filename=(str(in_dir)+'/'+str(coa_group)+'_2.fastq')

                    # Define Snakemake output files
                    output_files+=(path+"/"+final_temp_dir+"/"+coa_group+"_DASTool_bins ")

                    # Define Snakemake input files
                    if not (os.path.exists(coa1_filename) and os.path.exists(coa2_filename)):
                            # merge all .fastq for coassembly
                        merge1Cmd='cat '+read1_files+' > '+coa1_filename+''
                        subprocess.check_call(merge1Cmd, shell=True)

                        merge2Cmd='cat '+read2_files+' > '+coa2_filename+''
                        subprocess.check_call(merge2Cmd, shell=True)

                    else:
                        pass


        return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/coassembly_binning/Snakefile')

    # Run snakemake
    log_file=open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Metagenomics-Coassembly starting")
    log_file.close()

    mtg_snk_Cmd = 'snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    log_file=open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Metagenomics-Coassembly has finished :)")
    log_file.close()


    # Keep temp dirs / remove all
    if args.keep: # If -k, True: keep
        pass
    else: # If not -k, keep only last dir
        exist=list()
        for file in out_files.split(" "):
            exist.append(os.path.isfile(file))

        if all(exist): # all output files exist
            rmCmd='cd '+path+' | grep -v '+final_temp_dir+' | xargs rm -rf && mv '+final_temp_dir+' MCB_Holoflow'
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
