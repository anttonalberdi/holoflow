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
    config = os.path.join(os.path.abspath(curr_dir),"workflows/metagenomics/coassembly_binning/config.yaml")
else:
    config=args.config_file

if not (args.log):
    log = os.path.join(path,"Holoflow_coassembly_metagenomics.log")
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

    if data['assembler'] == "spades":
        merging=True
    else:
        merging=False


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
        # Paste desired output file names from input.txt
        read = 0
        group = 'empty'
        read1_files=''
        read2_files=''
        output_files=''
        final_temp_dir="MCB_04-BinMerging"

        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))

        for dir in lines:

            if not (dir.startswith('#')):
                dir = dir.strip('\n').split(' ') # Create a list of each line

                read+=1 # every sample will have two reads, keep the name of the file but change the read

                # Depending on spades or megahit, create a big file where all .fastq merged or concatenate by ,
                input_groupdir=str(dir[1])      # current input file path and name

                # Snakemake input files
                coa1_filename=(str(in_dir)+'/'+str(dir[0])+'_1.fastq')
                coa2_filename=(str(in_dir)+'/'+str(dir[0])+'_2.fastq')

                if merging: # spades is selected assembler
                        # write output files and finish group input
                    if group == 'empty': # will only happen on the first round - first coassembly group
                        group=dir[0]

                    elif ((not (group == dir[1])) or (line == last_line)): # when the group changes, define output files for previous group and finish input
                        #same as last output in Snakefile
                        output_files+=(path+"/"+final_temp_dir+"/"+group+"_DASTool_bins ")

                        # merge all .fastq for coassembly with spades
                        merge1Cmd='cd '+input_groupdir+' && cat *_1.fastq > '+coa1_filename+''
                        subprocess.check_call(merge1Cmd, shell=True)

                        merge2Cmd='cd '+input_groupdir+' && cat *_2.fastq > '+coa2_filename+''
                        subprocess.check_call(merge2Cmd, shell=True)

                        group=dir[0] # define new group in case first condition



                if not merging:   #megahit is the selected assembler, all files in string , separated

                        # write output files and finish group input
                    if group == 'empty': # will only happen on the first round - first coassembly group
                        group=dir[0]

                    elif ((not (group == dir[1])) or (line == last_line)): # when the group changes, define output files for previous group and finish input
                        #same as last output in Snakefile
                        output_files+=(path+"/"+final_temp_dir+"/"+group+"_DASTool_bins ")

                        # the .fastq files for megahit will contain a list of input files , separated instead of the read content
                        find1Cmd='find '+input_groupdir+'/*_1.fastq | tr "\n" "," > '+coa1_filename+''
                        subprocess.check_call(merge1Cmd, shell=True)

                        find2Cmd='find '+input_groupdir+'/*_2.fastq | tr "\n" "," > '+coa2_filename+''
                        subprocess.check_call(merge2Cmd, shell=True)

                        group=dir[0] # define new group in case first condition

        return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/coassembly_binning/Snakefile')

    # Run snakemake
    mtg_snk_Cmd = 'module unload gcc && module load tools anaconda3/4.4.0 && snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    print("Have a nice run!\n\t\tHOLOFOW Metagenomics-Coassembly starting")



###########################
#### Workflows running
###########################
# 2    # Metagenomics workflow
run_metagenomics(in_f, path, config, cores)
