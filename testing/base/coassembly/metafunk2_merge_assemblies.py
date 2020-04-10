#!/usr/bin/env python
# -*- coding: utf-8

"""The script for merging assemblies"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip
import signal

def merge_assemblies(projectname,projectpath,threads,memory,logfilepath):
    #Create merged and merged/assembly directories if do not exist
    merged_dir = "merged"
    merged_abs = os.path.join(projectpath, merged_dir)
    if not os.path.exists(merged_abs):
        os.makedirs(merged_abs)
    assembly_dir = "reassembly"
    assembly_abs = os.path.join(merged_abs, assembly_dir)
    if not os.path.exists(assembly_abs):
        os.makedirs(assembly_abs)

    assembliespath = os.path.join(projectpath,'*.assembly', 'contigs.fasta')
    mergedassembliespath = os.path.join(assembly_abs, 'assemblies.fna')
    mergedassembliesbase = os.path.join(assembly_abs, 'assemblies')
    nrassembliespath = os.path.join(assembly_abs, 'assemblies.nr.fna')
    afgassembliespath = os.path.join(assembly_abs, 'assemblies.afg')

    #Concatenate assemblies
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Concatenating assemblies \r\n".format(current_time))
    logfile.close()
    concCmd = 'cat '+assembliespath+' > '+mergedassembliespath+''
    subprocess.check_call(concCmd, shell=True)

    #Decomposed minimus2 pipeline (#https://github.com/nathanhaigh/amos/blob/master/src/Pipeline/minimus2.acf)
    #Consider using the alternative minimus2-blat - might be faster

    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running minimus2 pipeline to merge assemblies \r\n".format(current_time))
    logfile.close()

    mergedassemblies_bnk = os.path.join(mergedassembliesbase + '.bnk')
    mergedassemblies_afg = os.path.join(mergedassembliesbase + '.afg')
    mergedassemblies_refseq = os.path.join(mergedassembliesbase + '.ref.seq')
    mergedassemblies_qryseq = os.path.join(mergedassembliesbase + '.qry.seq')
    mergedassemblies_delta = os.path.join(mergedassembliesbase + '.delta')
    mergedassemblies_coords = os.path.join(mergedassembliesbase + '.coords')
    mergedassemblies_ovl = os.path.join(mergedassembliesbase + '.ovl')
    mergedassemblies_OVL = os.path.join(mergedassembliesbase + '.OVL')
    mergedassemblies_contig = os.path.join(mergedassembliesbase + '.contig')
    mergedassemblies_reassembly  = os.path.join(merged_abs, 'reassembly.fna')

    #Load software
    loadSoftware = 'module load perl/5.20.2 ncbi-blast/2.6.0+ cd-hit/4.8.1 MUMmer/3.23 kentUtils/350 amos/20121115 &&'

    #Modify merged assembly to afg format
    toamosCmd = ''+loadSoftware+' toAmos -s '+mergedassembliespath+' -o '+afgassembliespath+''
    subprocess.check_call(toamosCmd, shell=True)

    #Remove path if does not exist
    rmbankCmd = 'rm -fr '+mergedassemblies_bnk+''
    subprocess.check_call(rmbankCmd, shell=True)

    #Create bank
    bankCmd = ''+loadSoftware+' bank-transact -c -z -b '+mergedassemblies_bnk+' -m '+mergedassemblies_afg+''
    subprocess.check_call(bankCmd, shell=True)

    #Dump1
    dump1Cmd = ''+loadSoftware+' dumpreads '+mergedassemblies_bnk+' -M 0 > '+mergedassemblies_refseq+''
    subprocess.check_call(dump1Cmd, shell=True)

    #Dump2
    dump2Cmd = ''+loadSoftware+' dumpreads '+mergedassemblies_bnk+' -m 0 > '+mergedassemblies_qryseq+''
    subprocess.check_call(dump2Cmd, shell=True)

    #Nucmer
    nucmerCmd = ''+loadSoftware+' nucmer -maxmatch -c 100 '+mergedassemblies_refseq+' '+mergedassemblies_qryseq+' -p '+mergedassembliesbase+''
    subprocess.check_call(nucmerCmd, shell=True)

    #Coords
    coordsCmd = ''+loadSoftware+' show-coords -H -c -l -o -r -I 95 '+mergedassemblies_delta+' | nucmerAnnotate | egrep "BEGIN|END|CONTAIN|IDENTITY" > '+mergedassemblies_coords+''
    subprocess.check_call(coordsCmd, shell=True)

    #ovl
    ovlCmd = ''+loadSoftware+' nucmer2ovl -ignore 20 -tab '+mergedassemblies_coords+' | sort2 > '+mergedassemblies_ovl+''
    subprocess.check_call(ovlCmd, shell=True)

    #OVL
    OVLCmd = ''+loadSoftware+' ovl2OVL '+mergedassemblies_ovl+' > '+mergedassemblies_OVL+''
    subprocess.check_call(OVLCmd, shell=True)

    #Transact
    transactCmd = ''+loadSoftware+' bank-transact -z -b '+mergedassemblies_bnk+' -m '+mergedassemblies_OVL+''
    subprocess.check_call(transactCmd, shell=True)

    #Tigger
    tiggerCmd = ''+loadSoftware+' tigger -b '+mergedassemblies_bnk+''
    subprocess.check_call(tiggerCmd, shell=True)

    #Consensus
    consensusCmd = ''+loadSoftware+' make-consensus -B -e 0.06 -b '+mergedassemblies_bnk+' -w 15'
    subprocess.check_call(consensusCmd, shell=True)

    #Contig
    contigCmd = ''+loadSoftware+' bank2contig '+mergedassemblies_bnk+' > '+mergedassemblies_contig+''
    subprocess.check_call(contigCmd, shell=True)

    #Fasta
    fastaCmd = ''+loadSoftware+' bank2fasta -b '+mergedassemblies_bnk+' > '+mergedassemblies_reassembly+''
    subprocess.check_call(fastaCmd, shell=True)
