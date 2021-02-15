#!/bin/bash

infile=$1

if [ "$infile" == "" ] ; then
    echo "Usage: prokkagff2gtf.sh <PROKKA gff file>"
    exit 0
fi

grep -v "#" $infile | grep "UniProtKB" | sed -e 's/.*UniProtKB:\(.*\);locus.*/\1/' | sed -e 's/\$/\n/g' > UNIPROT

grep -v "#" $infile | grep "UniProtKB" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7 | sed -e 's/\$/\n/g' > PROKKA


paste PROKKA UNIPROT | awk -v OFS='\t' '{print $1,"PROKKA","CDS",$2,$3,".",$4,".","gene_id " $5}' && rm UNIPROT PROKKA
