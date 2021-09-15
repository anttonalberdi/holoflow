################################################################################
################################################################################
#
#Simple BASH script to automate holoflow project folder creation & CO job scripts
#Raphael Eisenhofer, September 2021
#
################################################################################
################################################################################

# Read in variables
echo "Enter directory name"
read directory

echo "Enter user name"
read user

echo "Welcome $user, please wait while I setup your holoflow project directory"

################################################################################

# Make directories
mkdir $directory
mkdir $directory/0_Configs
mkdir $directory/1_Scripts
mkdir $directory/2_References
mkdir $directory/3_Reports
mkdir $directory/logs
mkdir $directory/workdir

################################################################################

# Setup the preparegenomes script
echo "
#Declare full path to the project directory (the .sh file will be stored here as well)
projectpath=/home/projects/ku-cbd/people/$user/$directory
#Declare full path to holoflow
holoflowpath=/home/projects/ku-cbd/people/$user/holoflow
#Run holoflow
python3 $"{"holoflowpath"}"/preparegenomes.py \
-f $"{"projectpath"}"/0_Configs/preparegenomes_input.txt \
-d $"{"projectpath"}"/workdir \
-l $"{"projectpath"}"/logs/preparegenomes_$directory.log \
-c $"{"projectpath"}"/0_Configs/preparegenomes_config.yaml \
-t 40
" > $directory/1_Scripts/preparegenomes.sh

################################################################################

# Setup the preprocessing script
echo "
#Declare full path to the project directory (the .sh file will be stored here as well)
projectpath=/home/projects/ku-cbd/people/$user/$directory
#Declare full path to holoflow
holoflowpath=/home/projects/ku-cbd/people/$user/holoflow
#Run holoflow
python3 $"{"holoflowpath"}"/preprocessing.py \
-f $"{"projectpath"}"/0_Configs/preprocess_input.txt \
-d $"{"projectpath"}"/workdir \
-l $"{"projectpath"}"/logs/preprocess_$directory_log.log \
-adapter1 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA' \
-adapter2 'AAGTCGGATCGTAGCCATGTCGTTC' \
-N $directory_preprocess \
-g $"{"projectpath"}"/workdir/PRG/$directory_Human.fna \
-t 40
" > $directory/1_Scripts/preprocess.sh

################################################################################

# Setup the qsub script to send preparegenomes.sh to computerome job queue
echo "
qsub -V -A ku-cbd -W group_list=ku-cbd \
-d '`'pwd'`' \
-e $"{"projectpath"}"/logs/preparegenomes_$directory_error_file \
-o $"{"projectpath"}"/logs/preparegenomes_$directory_out_file.out \
-l nodes=1:ppn=40,mem=100gb,walltime=00:12:00:00 \
-N $directory_prepare_genomes \
$"{"projectpath"}"/1_Scripts/preparegenomes.sh
" > $directory/1_Scripts/submit_preparegenomes.sh

################################################################################

# Setup the qsub script to send preprocess.sh to computerome job queue
echo "
qsub -V -A ku-cbd -W group_list=ku-cbd \
-d '`'pwd'`' \
-e $"{"projectpath"}"/logs/preprocess_$directory_job_error_file \
-o $"{"projectpath"}"/logs/preprocess_$directory_job_out_file.out \
-l nodes=1:ppn=40,mem=100gb,walltime=00:48:00:00 \
-N $directory_preprocess \
$"{"projectpath"}"/1_Scripts/preprocess.sh
" > $directory/1_Scripts/submit_preprocess.sh

################################################################################
