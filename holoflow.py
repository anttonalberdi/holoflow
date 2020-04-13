import argparse
import subprocess
import os
import sys

###########################
#Argument parsing
###########################
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="input", dest="input", required=True)
parser.add_argument('-d', help="project directory path", dest="path", required=True)
parser.add_argument('-w', help="chosen workflow", dest="workflow", required=True)
parser.add_argument('-c', help="config file", dest="config", required=True)
args = parser.parse_args()

input=args.input
path=args.path
workflow=args.workflow
config=args.config



###########################
## Functions for input_output definition
###########################

def in_out_preprocessing(path,input):
    # Create "00-RawData/" directory if not exists
    in_dir = os.path.join(path,"00-InputData")
    if not os.path.exists(in_dir):
        os.makedirs(in_dir)

        with open(input,'r') as in_file:
            # Paste desired output file names from input.txt
            read = 0
            output_files=''

            lines = in_file.readlines()
            for file in lines:

                if not (file.startswith('#')):
                    file = file.strip('\n').split(' ')

                    read+=1
                    output_files+=(path+"/"+file[3]+"/"+file[0]+"_"+str(read)+".fastq ")

                    #Move files to new dir "00-RawData/" and change file names for 1st column in input.txt
                    filename=file[2]
                    copyfilesCmd='cp '+filename+' '+in_dir+'/'+file[0]+'_'+str(read)+'.fastq.gz'
                    subprocess.check_call(copyfilesCmd, shell=True)

                    if read == 2:
                        read=0
                        # Add stats output only once per sample
                        output_files+=(path+"/"+file[3]+"/"+file[0]+".stats ")

            return output_files

###########################
#### Snakemake pipeline run
###########################
load_modulesCmd='module unload gcc/5.1.0 && module load anaconda3/4.4.0'
subprocess.check_call(load_modulesCmd, shell=True)



###########################
#### Workflows
###########################

# 1    # Preprocessing workflow
if workflow == "preprocessing":

    # Define output names
    out_files = in_out_preprocessing(path,input)
    print(out_files)

    # Create preprocessing.sh for later job submission

    with open('./workflows/preprocessing/preprocessing.sh','w+') as sh:
        curr_dir = os.getcwd()
        path_snkf = os.path.join(curr_dir,'workflows/preprocessing/Snakefile')
        prep_snk = 'snakemake -s '+path_snkf+' '+out_files+' --configfile '+config+''
        sh.write(prep_snk)


    # Submit snakemake job
    preprocessingCmd = 'qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e '+path+'/Holo-preprocessing.err -o '+path+'/Holo-preprocessing.out -l nodes=1:ppn=40,mem=180gb,walltime=5:00:00:00 -N Holoflow-preprocessing ./workflows/preprocessing/preprocessing.sh'
    subprocess.check_call(preprocessingCmd, shell=True)
    print("Preprocessing with Holoflow was successfully submited")


# 2    # Metagenomics workflow

# if workflow == "metagenomics":
#
#     prep = input("Input files for holoflow/metagenomics are fastq. Is your data preprocessed? [y/n]")
#
#     if prep == 'n':
#         prep2 = input("Would you like to process it before running holoflow/metagenomics with holoflow/preprocessing? [y/n]")
#
#         if prep2 == 'n':
#             print("You should come back when your data is preprocessed. See you soon :)")
#         if prep2 == 'y':
#             snakemakeCmd = 'xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` '+path+'/snakemake.log -l nodes=1:ppn=28,mem=100gb,walltime=0:06:00:00 -N holoflow_metagenomics -de snakemake -s workflows/metagenomics/prep_and_metagenomics/Snakefile '+output_files+' --config '+config+''
#             subprocess.check_call(snakemakeCmd, shell=True)
#
#     if prep == 'y':
#         print("Great! Have a nice run!\n\t\tHOLOFOW Metagenomics starting")
#         snakemakeCmd = 'xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` '+path+'/snakemake.log -l nodes=1:ppn=28,mem=100gb,walltime=0:06:00:00 -N holoflow_metagenomics -de snakemake -s workflows/metagenomics/Snakefile '+output_files+' --config '+config+''
#         subprocess.check_call(snakemakeCmd, shell=True)


    # Genomics workflow
