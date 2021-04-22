in_data_drep=$1
in_data_checkm=$2
summary_table_tmp=$3
summary_table=$4



touch $summary_table_tmp
while read line; do
grep $line $in_data_drep | cut -d’,' -f1,2,3,5,6 >> $summary_table_tmp
done < <(cut -d’,' -f1 $in_data_checkm)
sort -t’,' -k2,2nr -k3,3n -k5,5nr $summary_table_tmp > $summary_table
rm $summary_table_tmp
#All MAGs
cat $summary_table | wc -l
#Near complete
awk -F ‘,’ ‘($2 > 98) && ($3 < 5) { print}’ $summary_table_tmp | wc -l
#High quality
awk -F ‘,’ ‘($2 > 90) && ($3 < 5) { print}’ $summary_table_tmp | wc -l
#Good quality
awk -F ‘,’ ‘($2 > 80) && ($3 < 10) { print}’ $summary_table_tmp | wc -l
