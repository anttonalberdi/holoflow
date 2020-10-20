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
parser.add_argument('-c', help="config file", dest="config_file", required=True)
parser.add_argument('-l', help="pipeline log file", dest="log", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
config=args.config_file
log=args.log
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
  data['logpath'] = str(log)
  dump = yaml.dump(data, config_file)



###########################
## Functions
###########################



    ###########################
    ###### PREPROCESSING FUNCTIONS

def in_out_preprocessing(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    # Define input directory and create it if not exists "00-InputData"
    in_dir = os.path.join(path,"PPR_00-InputData")
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Generate desired output file names from input.txt
        read = 0
        output_files=''
        final_temp_dir="PPR_03-MappedToReference"

        lines = in_file.readlines() # Read input.txt lines
        for file in lines:

            if not (file.startswith('#')):
                file = file.strip('\n').split(' ') # Create a list of each line

                read+=1     # every sample will have two reads, keep the name of the file but change the read
                # Add an output file based on input.txt info to a list for Snakemake command
                output_files+=(path+"/"+final_temp_dir+"/"+file[0]+"_"+str(read)+".fastq ")

                # Move files to new dir "00-InputData" and change file names for 1st column in input.txt
                #   if the current input file names do not match the designed ones in input.txt
                filename=file[2]       # current input file path and name
                desired_filename='"'+in_dir+'/'+file[0]+'_'+str(read)+'.fastq"'  # desired input file path and name specified in input.txt

                if not ((filename == desired_filename) and (os.path.exists(str(desired_filename)))):
                    if filename.endswith('.gz'):    # uncompress input file if necessary
                        uncompressCmd='gunzip -c '+filename+' > '+desired_filename+''
                        subprocess.check_call(uncompressCmd, shell=True)
                    else:                           # else just move the input file to "00-InputData" with the new name
                        copyfilesCmd='cp '+filename+' '+desired_filename+''
                        subprocess.check_call(copyfilesCmd, shell=True)


                if read == 2:
                    read=0  # two read files for one sample finished, new sample

                    # Add stats and bam output files only once per sample
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+".stats ")
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+"_ref.bam ")

        return output_files



################## NOT YET ###################

def prepare_threads(path,config):
    """Set a maximum number of used threads by AdapterRemoval during
    the quality filtering step based on the size and number of the
    input files"""

    # get input files average size:
    in_dir = os.path.join(path,"PPR_00-InputData")
    count_file_size=0
    number_file=0

    for file in os.listdir(in_dir):
        print(file)
        full_file=(''+in_dir+'/'+file+'')
        print(full_file)
        count_file_size+=os.path.getsize(os.path.abspath(full_file))
        number_file+=1

    # get average file size
    average_file_size = count_file_size/number_file
    number_file = number_file/2 # We count samples

    # depending on the average file size and number of files,
    # change number of threads for AdapterRemoval in config
    yaml = ruamel.yaml.YAML()
    yaml.explicit_start = True
    with open(str(config), 'r') as config_file:
        data = yaml.load(config_file)


    # If files smaller then 800MG, then it does not matter the num of samples w/4 threads
        # if files between 800MG and 1G then max 24 samples for 4 threads
        # if files between 1G and 2,5G then max 12 samples for 4 threads
        # if files between 2,5G and 5G then max 6 samples for 4 threads
        # if files between 5G and 10G then max 6 samples for 4 threads
    if (average_file_size < 800000000) or ((800000001 <= average_file_size <= 1000000000) and (number_file <= 24)) or ((1000000001 <= average_file_size <= 2500000000) and (number_file <= 12)) or ((2500000001 <= average_file_size <= 5000000000) and (number_file <= 6)) or ((5000000001 <= average_file_size <= 12000000000) and (number_file <= 3)):

        with open(str(config), 'w') as config_file:
            data['AdapterRemoval_threads'] = 4
            dump = yaml.dump(data, config_file)

    # Same corollary
    if ((800000001 <= average_file_size <= 1000000000) and (number_file > 24)) or ((1000000001 <= average_file_size <= 2500000000) and (12 < number_file <= 24)) or ((2500000001 <= average_file_size <= 5000000000) and (6 < number_file <= 12)) or ((5000000001 <= average_file_size <= 12000000000) and (3 < number_file <= 6)):

        with open(str(config), 'w') as config_file:
            data['AdapterRemoval_threads'] = 8
            dump = yaml.dump(data, config_file)

    # Same corollary
    if ((1000000001 <= average_file_size <= 2500000000) and (number_file > 24)) or ((2500000001 <= average_file_size <= 5000000000) and (12 < number_file <= 20)) or ((5000000001 <= average_file_size <= 12000000000) and (6 < number_file <= 10)):

        with open(str(config), 'w') as config_file:
            data['AdapterRemoval_threads'] = 14
            dump = yaml.dump(data, config_file)

    # Same corollary
    if ((2500000001 <= average_file_size <= 5000000000) and (number_file > 20)) or ((5000000001 <= average_file_size <= 10000000000) and (number_file > 10)):

        with open(str(log), 'w') as log_file:
            log_file.write("Your files are too big to be processed all together.\nIf these are average 12G, process maximum 10 files at a time.\nIf these are average 5G, process maximum 20 files at a time.")




def run_preprocessing(in_f, path, config, cores):
    """Run snakemake on shell"""

    # Define output names
    out_files = in_out_preprocessing(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/preprocessing/Snakefile')

    # get threads for AdapterRemoval
    prepare_threads(path,config)

    # Run snakemake
    prep_snk_Cmd = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.check_call(prep_snk_Cmd, shell=True)
    print("Have a nice run!\n\t\tHOLOFOW Preprocessing starting")

###########################
#### Snakemake pipeline run - load required modules
###########################
load_modulesCmd='module unload gcc && module load tools anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)



###########################
#### Workflows running
###########################


# 1    # Preprocessing workflow

run_preprocessing(in_f, path, config, cores)
