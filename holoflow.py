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
parser.add_argument('-w', help="chosen workflow", dest="workflow", required=True)
parser.add_argument('-c', help="config file", dest="config_file", required=True)
parser.add_argument('-t', help="threads", dest="threads", required=True)
args = parser.parse_args()

in_f=args.input_txt
path=args.work_dir
workflow=args.workflow
config=args.config_file
cores=args.threads



###########################
## Functions
###########################

    ###########################
    ###### PREPROCESSING FUNCTIONS

def in_out_preprocessing(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    # Create "00-RawData/" directory if not exists
    in_dir = os.path.join(path,"00-InputData")
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Paste desired output file names from input.txt
        read = 0
        output_files=''
        final_temp_dir="PPR04-MappedToHuman"

        lines = in_file.readlines()
        for file in lines:

            if not (file.startswith('#')):
                file = file.strip('\n').split(' ')

                read+=1
                output_files+=(path+"/"+final_temp_dir+"/"+file[0]+"_"+str(read)+".fastq ")

                #Move files to new dir "00-InputData" and change file names for 1st column in input.txt
                filename=file[2]
                desired_filename='"'+in_dir+'/'+file[0]+'_'+str(read)+'.fastq.gz"'
                if not (filename == desired_filename):
                    copyfilesCmd='cp '+filename+' '+desired_filename+''
                    subprocess.check_call(copyfilesCmd, shell=True)

                if read == 2:
                    read=0
                    # Add stats output only once per sample
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+".stats ")

        return output_files


def run_preprocessing(in_f, path, config, cores):
    """Create preprocessing.sh file and run snakemake on shell"""
    # Define output names
    out_files = in_out_preprocessing(path,in_f)

    # Create preprocessing.sh for later job submission
    with open('./workflows/preprocessing/preprocessing.sh','w+') as sh:
        curr_dir = os.getcwd()
        path_snkf = os.path.join(curr_dir,'workflows/preprocessing/Snakefile')
        prep_snk = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+' --cores '+cores+''
        sh.write(prep_snk)


    # Submit snakemake job
    preprocessingCmd = 'qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e '+path+'/Holo-preprocessing.err -o '+path+'/Holo-preprocessing.out -l nodes=1:ppn=40,mem=180gb,walltime=5:00:00:00 -N Holoflow-preprocessing ./workflows/preprocessing/preprocessing.sh'
    subprocess.check_call(preprocessingCmd, shell=True)
    print("Preprocessing with Holoflow was successfully submited")



    ###########################
    ###### METAGENOMICS FUNCTIONS

def in_out_metagenomics(path,in_f):
    """Generate output names files from input.txt. Rename and move
    input files where snakemake expects to find them if necessary."""
    in_dir = os.path.join(path,"04-MappedToHuman")
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

    with open(in_f,'r') as in_file:
        # Paste desired output file names from input.txt
        read = 0
        output_files=''
        final_temp_dir="MIA03-Binning"

        lines = in_file.readlines()
        for file in lines:

            if not (file.startswith('#')):
                file = file.strip('\n').split(' ')

                read+=1

                output_files+=(path+"/"+final_temp_dir+"/"+file[0]+"_dastool")


                #Move files to input dir "04-MappedToHuman/" and change file names for column 1 in input.txt
                filename=file[2]
                desired_filename='"'+in_dir+'/'+file[0]+'_'+str(read)+'.fastq.gz"'

                if not (filename == desired_filename):
                    copyfilesCmd='cp '+filename+' '+desired_filename+''
                    subprocess.check_call(copyfilesCmd, shell=True)

                if read == 2:
                    read=0
                    # Add stats output only once per sample
                    output_files+=(path+"/"+final_temp_dir+"/"+file[0]+".stats ")

        return output_files



def run_metagenomics(in_f, path, config, cores):
    """Create metagenomics.sh file and run snakemake on shell"""

    # Define output names
    out_files = in_out_metagenomics(path,in_f)


    # # Create preprocessing.sh for later job submission
    with open('./workflows/metagenomics/metagenomics.sh','w+') as sh:
        curr_dir = os.getcwd()
        path_snkf = os.path.join(curr_dir,'workflows/metagenomics/Snakefile')
        prep_snk = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+' --cores '+cores+''
        sh.write(prep_snk)


    metagenomicsCmd = 'qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e '+path+'/Holo-metagenomics.err -o '+path+'/Holo-metagenomics.out -l nodes=1:ppn=40,mem=180gb,walltime=5:00:00:00 -N Holoflow-metagenomics ./workflows/metagenomics/metagenomics.sh'
    subprocess.check_call(metagenomicsCmd, shell=True)
    print("Great! Have a nice run!\n\t\tHOLOFOW Metagenomics starting")



    ###########################
    ###### PREPROCESSING AND METAGENOMICS FUNCTIONS

# def run_prepandmet(prepin_f, metin_f, path, prepconfig, metconfig, cores):
#     """Run both preprocessing and metagenomics Snakefiles on shell"""
#
#     # Define output names
#     out_files = in_out_metagenomics(path,in_f)
#


###########################
#### Snakemake pipeline run - load required modules
###########################
load_modulesCmd='module unload gcc/5.1.0 && module load anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)



###########################
#### Workflows
###########################

# 1    # Preprocessing workflow
if workflow == "preprocessing":
    run_preprocessing(in_f, path, config, cores)


# 2    # Metagenomics workflow

if workflow == "metagenomics":

    prepdata = input("Is your data preprocessed into fastq files? [y/n]")

    if prepdata == 'n':
        prepdata2 = input("Would you like to process it before running holoflow/metagenomics with holoflow/preprocessing? [y/n]")

        if prepdata2 == 'n':
            print("You should come back when your data is preprocessed. See you soon :)")

        # if prepdata2 == 'y':    # It would be much easier to concatenate Snakefiles and new functions - DO IT
        #     prep_in_f = input("Could you please state the path for the preprocessing input file? - No quoting needed\n")
        #     prep_config = input("Could you please state the path for the preprocessing config file? - No quoting needed\n")
        #     run_preprocessing(prep_in_f, path, prep_config, cores)
        #
        #     prep_out_dir = os.path.join(path,"04-MappedToHuman")
        #     if os.path.exists(prep_out_dir):
        #         run_metagenomics(in_f, path, config, cores)

    if prepdata == 'y':
        run_metagenomics(in_f, path, config, cores)




# 3    # Genomics workflow
