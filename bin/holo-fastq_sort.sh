# Sort fastq files
fastq1="/home/projects/ku-cbd/people/nurher/coassembly_test_BATS/PPR_03-MappedToReference/Bats_coa_groupB_1.fastq"
fastq2="/home/projects/ku-cbd/people/nurher/coassembly_test_BATS/PPR_03-MappedToReference/Bats_coa_groupB_2.fastq"
sortedfq1="/home/projects/ku-cbd/people/nurher/coassembly_test_BATS/PPR_03-MappedToReference/Sorted_groupB_1.fastq"
sortedfq2="/home/projects/ku-cbd/people/nurher/coassembly_test_BATS/PPR_03-MappedToReference/Sorted_groupB_2.fastq"


cat ${fastq1} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sortedfq1}
cat ${fastq2} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sortedfq2}
