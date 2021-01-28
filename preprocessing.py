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
parser.add_argument('-g', help="reference genome", dest="ref", required=False)
parser.add_argument('-k', help="keep tmp directories", dest="keep", action='store_true')
parser.add_argument('-l', help="pipeline log file", dest="log", required=False)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
ref=args.ref
cores=args.threads

    # retrieve current directory
file = os.path.dirname(sys.argv[0])
curr_dir = os.path.abspath(file)


if not (args.config_file):
    config = os.path.join(os.path.abspath(curr_dir),"workflows/preprocessing/config.yaml")
else:
    config=args.config_file

if not (args.log):
    log = os.path.join(path,"Holoflow_preprocessing.log")
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
    data['refgenomes'] = str(ref)
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
        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines)))

        # Define variables
        output_files=''
        final_temp_dir="PPR_03-MappedToReference"

        for line in lines:
            ### Skip line if starts with # (comment line)
            if not (line.startswith('#')):

                line = line.strip('\n').split(' ') # Create a list of each line
                sample_name=line[0]
                in_for=line[1]
                in_rev=line[2]

                # Define output files based on input.txt
                output_files+=path+'/'+final_temp_dir+'/'+sample_name+'_1.fastq '
                output_files+=path+'/'+final_temp_dir+'/'+sample_name+'_2.fastq '


                # Define input file
                in1=in_dir+'/'+sample_name+'_1.fastq.tmp'
                # Check if input files already in desired dir
                if os.path.isfile(in1):
                    pass
                else:
                    #If the file is not in the working directory, transfer it
                    if os.path.isfile(in_for) and not (os.path.isfile(in1)):
                        if in_for.endswith('.gz'):
                            read1Cmd = 'ln -s '+in_for+' '+in1+'.gz && gunzip -c '+in1+'.gz > '+in1+''
                            subprocess.Popen(read1Cmd, shell=True).wait()
                        else:
                            read1Cmd = 'ln -s '+in_for+' '+in1+''
                            subprocess.Popen(read1Cmd, shell=True).wait()


                # Define input file
                in2=in_dir+'/'+sample_name+'_2.fastq.tmp'
                # Check if input files already in desired dir
                if os.path.isfile(in2):
                    pass
                else:
                    #If the file is not in the working directory, transfer it
                    if os.path.isfile(in_rev) and not (os.path.isfile(in2)):
                        if in_for.endswith('.gz'):
                            read2Cmd = 'ln -s '+in_rev+' '+in2+'.gz && gunzip -c '+in2+'.gz > '+in2+''
                            subprocess.Popen(read2Cmd, shell=True).wait()
                        else:
                            read2Cmd = 'ln -s '+in_rev+' '+in2+''
                            subprocess.Popen(read2Cmd, shell=True).wait()


                # Add stats and bam output files only once per sample
                output_files+=(path+"/"+final_temp_dir+"/"+sample_name+".stats ")
                output_files+=(path+"/"+final_temp_dir+"/"+sample_name+"_ref.bam ")

        return output_files



def run_preprocessing(in_f, path, config, cores):
    """Run snakemake on shell, wait for it to finish.
    Given flag, decide whether keep only last directory."""

    # Define output names
    out_files = in_out_preprocessing(path,in_f)
    curr_dir = os.path.dirname(sys.argv[0])
    holopath = os.path.abspath(curr_dir)
    path_snkf = os.path.join(holopath,'workflows/preprocessing/Snakefile')

    # Run snakemake
    log_file = open(str(log),'w+')
    log_file.write("Have a nice run!\n\t\tHOLOFOW Preprocessing starting")
    log_file.close()

    prep_snk_Cmd = 'module load tools anaconda3/4.4.0 && snakemake -s '+path_snkf+' -k '+out_files+' --configfile '+config+' --cores '+cores+''
    subprocess.Popen(prep_snk_Cmd, shell=True).wait()

    log_file = open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Preprocessing has finished :)")
    log_file.close()

    # Keep temp dirs / remove all
    if args.keep: # If -k, True: keep
        pass
    else: # If not -k, keep only last dir
        exist=list()
        for file in out_files.split(" "):
            exist.append(os.path.isfile(file))

        if all(exist): # all output files exist
            rmCmd='cd '+path+' | grep -v '+final_temp_dir+' | xargs rm -rf && mv '+final_temp_dir+' PPR_Holoflow'
            subprocess.Popen(rmCmd,shell=True).wait()

        else:   # all expected output files don't exist: keep tmp dirs
            log_file = open(str(log),'a+')
            log_file.write("Looks like something went wrong...\n\t\t The temporal directories have been kept, you should have a look...")
            log_file.close()




###########################
#### Workflows running
###########################


# 1    # Preprocessing workflow
run_preprocessing(in_f, path, config, cores)
