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
    config = os.path.join(os.path.abspath(curr_dir),"workflows/metagenomics/individual_binning/config.yaml")
else:
    config=args.config_file

if not (args.log):
    log = os.path.join(path,"Holoflow_individualA_metagenomics.log")
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


###########################
## Functions
###########################

    ###########################
    ###### METAGENOMICS FUNCTIONS

def in_out_metagenomics(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    in_dir = os.path.join(path,"PPR_03-MappedToReference")

    if os.path.exists(in_dir):
        rmdirCmd='cd '+in_dir+'/.. && rm -rf '+in_dir+' && mkdir '+in_dir+''
        subprocess.check_call(rmdirCmd,shell=True)

    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Define variables
        output_files=''
        final_temp_dir="MIB_04-BinMerging"
        all_lines = in_file.readlines() # Read input.txt lines

        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))

        for line in lines:
            ### Skip line if starts with # (comment line)
            if not (line.startswith('#')):

                line = line.strip('\n').split(' ') # Create a list of each line
                sample_name=line[0]
                in_for=line[1]
                in_rev=line[2]


                # Define input file
                in1=in_dir+'/'+sample_name+'_1.fastq'
                # Check if input files already in desired dir
                if os.path.isfile(in1):
                    pass
                else:
                    #If the file is not in the working directory, transfer it
                    if os.path.isfile(in_for):
                        if in_for.endswith('.gz'):
                            read1Cmd = 'gunzip -c '+in_for+' > '+in1+''
                            subprocess.Popen(read1Cmd, shell=True).wait()
                        else:
                            read1Cmd = 'cp '+in_for+' '+in1+''
                            subprocess.Popen(read1Cmd, shell=True).wait()


                # Define input file
                in2=in_dir+'/'+sample_name+'_2.fastq'
                # Check if input files already in desired dir
                if os.path.isfile(in2):
                    pass
                else:
                    #If the file is not in the working directory, transfer it
                    if os.path.isfile(in_rev):
                        if in_for.endswith('.gz'):
                            read2Cmd = 'gunzip -c '+in_rev+' > '+in2+''
                            subprocess.Popen(read2Cmd, shell=True).wait()
                        else:
                            read2Cmd = 'cp '+in_rev+' '+in2+''
                            subprocess.Popen(read2Cmd, shell=True).wait()


                output_files+=(path+"/"+final_temp_dir+"/"+sample_name+"_DASTool_bins ")

        return output_files




def run_metagenomics(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/metagenomics/individual_binning/Snakefile')

    # Run snakemake
    mtg_snk_Cmd = 'module unload gcc && module load tools anaconda3/4.4.0 && snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(mtg_snk_Cmd, shell=True)

    print("Have a nice run!\n\t\tHOLOFOW Metagenomics-IndividualBinning starting")



###########################
#### Workflows running
###########################
# 2    # Metagenomics workflow
run_metagenomics(in_f, path, config, cores)
