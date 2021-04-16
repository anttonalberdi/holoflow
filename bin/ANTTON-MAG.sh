touch ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.tmp.csv
while read line; do
grep $line ${workdir}/MDR_01-BinDereplication/${batch}/data_tables/genomeInformation.csv | cut -d’,' -f1,2,3,5,6 >> ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.tmp.csv
done < <(cut -d’,' -f1 MDR_01-BinDereplication/${batch}/data_tables/Widb.csv)
sort -t’,' -k2,2nr -k3,3n -k5,5nr ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.tmp.csv > ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.csv
rm ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.tmp.csv
#All MAGs
cat ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.csv | wc -l
#Near complete
awk -F ‘,’ ‘($2 > 98) && ($3 < 5) { print}’ ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.csv | wc -l
#High quality
awk -F ‘,’ ‘($2 > 90) && ($3 < 5) { print}’ ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.csv | wc -l
#Good quality
awk -F ‘,’ ‘($2 > 80) && ($3 < 10) { print}’ ${workdir}/MDR_01-BinDereplication/${batch}/derep_bins_Info.csv | wc -l
