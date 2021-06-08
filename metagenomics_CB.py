import argparse
import subprocess
import os
import re
import glob
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
parser.add_argument('-W', help="threads", dest="REWRITE", action='store_true')
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
    cpconfigCmd= 'cp '+curr_dir+'/workflows/metagenomics/coassembly_binning/config.yaml '+path+'/'+current_time+'_config.yaml'
    subprocess.Popen(cpconfigCmd,shell=True).wait()


    config = path+'/'+current_time+'_config.yaml'
else:
    config=args.config_file

# If the user does not specify a log file, provide default path
if not (args.log):
    log = os.path.join(path,"Holoflow_coassembly_metagenomics.log")
else:
    log=args.log


    # Load dependencies
loaddepCmd='module unload gcc && module load tools anaconda3/4.4.0'
subprocess.Popen(loaddepCmd,shell=True).wait()


    #Append variables to .yaml config file for Snakefile calling standalone files
import ruamel.yaml
yaml = ruamel.yaml.YAML() # create yaml obj
yaml.explicit_start = True
with open(str(config), 'r') as config_file:
    data = yaml.load(config_file)# get data found now in config - as dictionary
    if data == None: # if config is empty, create dictionary
        data = {}

with open(str(config), 'w') as config_file:
    data['threads'] = str(cores)
    data['holopath'] = str(curr_dir)
    data['logpath'] = str(log)
    dump = yaml.dump(data, config_file) # load updated dictionary to config file


###########################
## Functions
###########################

    ###########################
    ###### METAGENOMICS FUNCTIONS

def in_out_metagenomics(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    in_dir = os.path.join(path,"PPR_03-MappedToReference")
    merged_in_dir = os.path.join(path,"MCB_00-MergedData")

    if not os.path.exists(in_dir): # create dir with all files to input to co-assembly
        os.makedirs(in_dir)
    else:
        pass

    # create dir for merged files (2 files containing data of all inputted files)
    if not os.path.exists(merged_in_dir):
        os.makedirs(merged_in_dir)
    else:
        pass

    with open(in_f,'r') as in_file:
        # Define necessary variables
        coa_group = False # coassembly group ID still not defined
        coa1_filename=''
        coa2_filename=''
        read1_files=''
        list_read1=list()
        read2_files=''
        list_read2=list()
        output_files=''
        final_temp_dir="MCB_04-BinMerging"

        all_lines = in_file.readlines() # Read input.txt lines
        # remove empty lines
        all_lines = map(lambda s: s.strip(), all_lines)
        lines = list(filter(None, list(all_lines))) # save input file content withput blank lines in "lines"
        last_line = lines[-1].split(' ') # last line of input file


    for line in lines:

        if not (line.startswith('#')):
            line = line.strip('\n').split(' ') # Create a list of each line
            sample=str(line[0])      # sample ID


            if (coa_group == line[1]) or not coa_group: # If sample in same coa group or not defined yet

                read1_files+=line[2]+' '
                read2_files+=line[3]+' '
                coa_group=line[1]


            if coa_group and not (coa_group == line[1]): # When the coa group is defined and changes, define output files for previous group and finish input

                # Fill in PPR_03 of uniformely renamed files
                input_dir = in_dir+'/'+coa_group
                if os.path.exists(input_dir):
                    if args.REWRITE: # If user wants to remove previous runs' data and run from scratch
                        rmCmd='rm -rf '+input_dir+''
                        subprocess.Popen(rmCmd,shell=True).wait()
                    else:
                        pass
                if not os.path.exists(input_dir): # if input directory does not exist
                    os.makedirs(input_dir)


                ###### Handle individual sample files before merging them
                list_read1=read1_files.strip().split(' ')
                list_read2=read2_files.strip().split(' ')

                for file1 in list_read1:
                    file=os.path.basename(file1)
                    # fastq inputted files to coassembly can have various nomenclatures
                    # _1.fastq, _1.fq, .1.fastq, .1.fq, etc.
                    #This command retrieves the file ID without format and for/rev number
                    sampleID=re.sub('(\.|_)[0-9]{1}\.f[aA-zZ]*\.?.*','',file)

                    # create a standardized directory with standardized IDs to coassemble
                    if file1.endswith('.gz'):
                        read1=input_dir+'/'+sampleID+'_1.fastq.gz'
                    else:
                        read1=input_dir+'/'+sampleID+'_1.fastq'

                    try:
                        cp1Cmd='ln -s '+file1+' '+read1+''  # If the file already existed, won't create link
                        subprocess.Popen(cp1Cmd, shell=True).wait()
                    except:
                        pass

                for file2 in list_read2:
                    file=os.path.basename(file2)
                    sampleID=re.sub('(\.|_)[0-9]{1}\.f[aA-zZ]*\.?.*','',file)

                    if file2.endswith('.gz'):
                        read2=in_dir+'/'+coa_group+'/'+sampleID+'_2.fastq.gz'
                    else:
                        read2=in_dir+'/'+coa_group+'/'+sampleID+'_2.fastq'

                    try:
                        cp2Cmd='ln -s '+file2+' '+read2+''  # If the file already existed, won't create link
                        subprocess.Popen(cp2Cmd, shell=True).wait()
                    except:
                        pass

            ###### Create coassembly merged files from all individual samples
                coa1_filename=(str(merged_in_dir)+'/'+str(coa_group)+'_1.fastq')
                coa2_filename=(str(merged_in_dir)+'/'+str(coa_group)+'_2.fastq')

                # if the forward read merged file exists, choose if rewrite or not
                if os.path.isfile(coa1_filename):
                    if args.REWRITE: # If user wants to remove previous runs' data and run from scratch
                        rmCmd='rm '+coa1_filename+' '+coa2_filename+''
                        subprocess.Popen(rmCmd,shell=True).wait()
                    else: #user wants to continue from rpevious run
                        pass

                if not os.path.isfile(coa1_filename):
                    files1 = glob.glob(in_dir+'/'+coa_group+'/*_1.fastq*')
                    for file1 in files1:
                        # Create a files called ".fastq", but actually fill them with a comma-separarted
                        # string of all the files that want to be considered for the coassembly
                        # MEGAHIT accepts this string as input, while MetaSpades will require the actual
                        #  merging of the files into 1 file: done in holo-assembly file -> only for SMALL coassemblies!
                        with open(coa1_filename,'a+') as coa1, open(coa2_filename,'a+') as coa2:
                            if file1 == files1[-1]:
                                coa1.write(file1.strip())

                                file2 = file1.strip().replace('1.fastq','2.fastq')
                                coa2.write(file2.strip())
                            else:
                                coa1.write(file1.strip()+',')

                                file2 = file1.strip().replace('1.fastq','2.fastq')
                                coa2.write(file2.strip()+',')

                # Define Snakemake output files
                output_files+=(path+"/"+final_temp_dir+"/"+coa_group+"_files ")

                # Define new coa group
                coa_group=line[1]
                read1_files=''
                read1_files+=line[2]+' '
                list_read1=list()
                read2_files=''
                read2_files+=line[3]+' '
                list_read2=list()


            if line == last_line: # in this case it is as if the coassembly group was changing, finish
                # Fill in PPR_03 of uniformely renamed files
                input_dir = in_dir+'/'+coa_group
                if os.path.exists(input_dir):
                    if args.REWRITE:
                        rmCmd='rm -rf '+input_dir+''
                        subprocess.Popen(rmCmd,shell=True).wait()
                    else:
                        pass
                if not os.path.exists(input_dir):
                    os.makedirs(input_dir)


                ###### Handle individual sample files
                list_read1=read1_files.strip().split(' ')
                list_read2=read2_files.strip().split(' ')

                for file1 in list_read1:
                    file=os.path.basename(file1)
                    sampleID=re.sub('(\.|_)[0-9]{1}\.f[aA-zZ]*\.?.*','',file)

                    if file1.endswith('.gz'):
                        read1=in_dir+'/'+coa_group+'/'+sampleID+'_1.fastq.gz'
                    else:
                        read1=in_dir+'/'+coa_group+'/'+sampleID+'_1.fastq'

                    try:
                        cp1Cmd='ln -s '+file1+' '+read1+''  # If the file already existed, won't create link
                        subprocess.Popen(cp1Cmd, shell=True).wait()
                    except:
                        pass

                for file2 in list_read2:
                    file=os.path.basename(file2)
                    sampleID=re.sub('(\.|_)[0-9]{1}\.f[aA-zZ]*\.?.*','',file)

                    if file2.endswith('.gz'):
                        read2=in_dir+'/'+coa_group+'/'+sampleID+'_2.fastq.gz'
                    else:
                        read2=in_dir+'/'+coa_group+'/'+sampleID+'_2.fastq'

                    try:
                        cp2Cmd='ln -s '+file2+' '+read2+''  # If the file already existed, won't create link
                        subprocess.Popen(cp2Cmd, shell=True).wait()
                    except:
                        pass

            ###### Create coassembly files data
                coa1_filename=(str(merged_in_dir)+'/'+str(coa_group)+'_1.fastq')
                coa2_filename=(str(merged_in_dir)+'/'+str(coa_group)+'_2.fastq')

                if os.path.isfile(coa1_filename):
                    if args.REWRITE:
                        rmCmd='rm '+coa1_filename+' '+coa2_filename+''
                        subprocess.Popen(rmCmd,shell=True).wait()
                    else:
                        pass

                if not os.path.isfile(coa1_filename):
                    files1 = glob.glob(in_dir+'/'+coa_group+'/*_1.fastq*')
                    for file1 in files1:
                        with open(coa1_filename,'a+') as coa1, open(coa2_filename,'a+') as coa2:
                            if file1 == files1[-1]:
                                coa1.write(file1.strip())

                                file2 = file1.strip().replace('1.fastq','2.fastq')
                                coa2.write(file2.strip())
                            else:
                                coa1.write(file1.strip()+',')

                                file2 = file1.strip().replace('1.fastq','2.fastq')
                                coa2.write(file2.strip()+',')

                # Define Snakemake output files
                output_files+=(path+"/"+final_temp_dir+"/"+coa_group+"_files ")

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
    subprocess.Popen(mtg_snk_Cmd, shell=True).wait()

    log_file=open(str(log),'a+')
    log_file.write("\n\t\tHOLOFOW Metagenomics-Coassembly has finished :)")
    log_file.close()


    # Keep temp dirs / remove all
    if args.keep: # If -k, True: keep
        pass
    else: # If not -k, keep only last dir
        exist=list()
        for file in out_files.split(' '):
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
