#!/bin/bash

# get data from stdin
env_file=$1
config_dbs=$2
assembly=$3
out_dir=$4
threads=$5
min_c_size=$6

module load anaconda3/4.4.0
conda env create -f $env_file -n DRAM
conda activate DRAM
DRAM-setup.py import_config --config_loc $config_dbs
DRAM.py annotate -i $assembly -o $out_dir --threads $threads --min_contig_size $min_c_size

# define vars for distill
in="${out_dir}/annotations.tsv"
out="${out_dir}/summary"
trna="${out_dir}/trnas.tsv"
rrna="${out_dir}/rrnas.tsv"

DRAM.py distill -i $in -o $out --trna_path $trna --rrna_path $rrna
