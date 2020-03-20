#!/usr/bin/env python
# -*- coding: utf-8

"""The script for contig binning from reassembly"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip
import glob
import shutil

def binning_merged(projectname,projectpath,threads,memory,logfilepath):
    binningdir = "binning"
    binningdir_abs = os.path.join(projectpath, 'merged',binningdir)
    if not os.path.exists(binningdir_abs):
        os.makedirs(binningdir_abs)

    #Declare input files
    reassemblypath = os.path.join(projectpath, 'merged', 'reassembly.fna')
    reassemblybampaths = os.path.join(projectpath, 'merged','reassembly_mapping','*.sorted.bam')

    #########################
    ######## Metabat ########
    #########################

    metabatdir = os.path.join(binningdir_abs, 'metabat')
    if not os.path.exists(metabatdir):
        os.makedirs(metabatdir)
    metabatdepthfile = os.path.join(metabatdir, 'depth.txt')
    metabatbinbase = os.path.join(metabatdir, 'metabat')

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating metabat depth file from the reads mapped to the reassembly \r\n".format(current_time))
    logfile.close()
    metabatdepthfileCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+metabatdepthfile+' '+reassemblybampaths+'' #added tools
    subprocess.check_call(metabatdepthfileCmd, shell=True)

    #Run metabat
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running metabat binning\r\n".format(current_time))
    logfile.close()
    metabatCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && metabat2 -i '+reassemblypath+' -a '+metabatdepthfile+' -o '+metabatbinbase+' -m 1500 -t '+threads+''
    subprocess.check_call(metabatCmd, shell=True)

    #Create contig to bin table
    bintablefile = os.path.join(binningdir_abs, 'bins_metabat.txt')
    bintable=open(bintablefile,"a+")
    metabatdir = os.path.join(metabatbinbase + '.*.fa')
    binlist = glob.glob(metabatdir)
    for bin in binlist:
        binname = os.path.splitext(os.path.basename(bin))[0]+''
        with open(bin, 'r') as binfile:
           for line in binfile:
                if line.startswith('>'):
                    contig = line.strip()
                    contig = contig.replace(">", "")
                    bintable.write("{0}\t{1}\r\n".format(contig,binname))
    bintable.close()

    #########################
    ######## Maxbin #########
    #########################

    maxbindir = os.path.join(binningdir_abs, 'maxbin')
    if not os.path.exists(maxbindir):
        os.makedirs(maxbindir)
    maxbindepthfile = os.path.join(maxbindir, 'depth.txt')
    maxbinbase = os.path.join(maxbindir, 'maxbin')

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating maxbin depth file from the reads mapped to the assembly \r\n".format(current_time))
    logfile.close()
    maxbindepthfileCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+maxbindepthfile+' --noIntraDepthVariance '+reassemblybampaths+''
    subprocess.check_call(maxbindepthfileCmd, shell=True)

    #Run maxbin
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running maxbin \r\n".format(current_time))
    logfile.close()
    maxbinCmd = 'module load perl/5.20.2 maxbin/2.2.7 fraggenescan/1.31 && run_MaxBin.pl -contig '+reassemblypath+' -abund '+maxbindepthfile+' -out '+maxbinbase+' -thread '+threads+''
    subprocess.check_call(maxbinCmd, shell=True)

    #Create contig to bin table
    bintablefile = os.path.join(binningdir_abs, 'bins_maxbin.txt')
    bintable=open(bintablefile,"a+")
    maxbindir = os.path.join(maxbinbase + '.*.fasta')
    binlist = glob.glob(maxbindir)
    for bin in binlist:
        binname = os.path.splitext(os.path.basename(bin))[0]+''
        with open(bin, 'r') as binfile:
           for line in binfile:
                if line.startswith('>'):
                    contig = line.strip()
                    contig = contig.replace(">", "")
                    bintable.write("{0}\t{1}\r\n".format(contig,binname))
    bintable.close()

def bin_refinement(projectname,projectpath,threads,memory,logfilepath):
    bincontig_tables = ",".join(glob.glob(os.path.join(projectpath,'merged/binning', 'bins_*.txt')))
    reassemblypath = os.path.join(projectpath, 'merged', 'reassembly.fna')
    dastoolpath = os.path.join(projectpath, 'merged','binning','dastool')
    if not os.path.exists(dastoolpath):
        os.makedirs(dastoolpath)
    dastoolbase = os.path.join(dastoolpath, 'dastool')

    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Refinning bins using DAS_Tool \r\n".format(current_time))
    logfile.close()

    #Refinement using DAS_Tool
    dastooldb = '/home/projects/ku-cbd/people/antalb/databases/dastool_db'
    dastoolDependencies = 'module load gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ prodigal/2.6.3 das_tool/1.1.1 diamond/0.9.24 usearch/11.0.667'
    dastoolCmd = ''+dastoolDependencies+' && DAS_Tool -i '+bincontig_tables+' -c '+reassemblypath+' -o '+dastoolbase+' -l maxbin,metabat --search_engine diamond -t '+threads+' --db_directory '+dastooldb+' --write_bins 1'
    subprocess.check_call(dastoolCmd, shell=True)

    #Refinement using Binning_refiner (problems with R dependencies)
    #module unload gcc gcc/5.1.0 && module load anaconda3/4.0.0 && Binning_refiner -i metafunk2_test2/merged/binning/refiner/ -p refined -plot

    #Move definitive bins to binning directory
    binsource = os.path.join(projectpath, 'merged','binning','dastool','dastool_DASTool_bins')
    bindestination = os.path.join(projectpath, 'merged','binning')
    binfiles = glob.glob(os.path.join(binsource,'*.fa'))
    for b in binfiles:
        shutil.move(b, bindestination)
