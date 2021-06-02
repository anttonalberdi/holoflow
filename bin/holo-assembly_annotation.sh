#!/bin/bash

# get data from stdin
env_file=$1       # environment file to create conda env
config_dbs=$2     # config file for DRAM to load and know where dbs were downloaded
assembly=$3
out_dir=$4
threads=$5
min_c_size=$6     # minimum config size

module load miniconda3/4.10.1
conda env create -f $env_file -n DRAM   # create conda environment
wait
#conda init bash
conda activate DRAM     # activate conda environment
DRAM-setup.py import_config --config_loc $config_dbs    # load config file for DRAM to know where to find dbs
DRAM.py annotate -i $assembly -o $out_dir --threads $threads --min_contig_size $min_c_size    # functional annotation

# define vars for distill
in="${out_dir}/annotations.tsv"
out="${out_dir}/summary"
trna="${out_dir}/trnas.tsv"
rrna="${out_dir}/rrnas.tsv"

DRAM.py distill -i $in -o $out --trna_path $trna --rrna_path $rrna  # extract most relevant info
wait
wait
conda deactivate    # deactivate conda environment
