in_data_drep=$1
in_data_checkm=$2
summary_table_tmp=$3
mag_table=$4
summary_table=$5


touch $summary_table_tmp
while read line; do
grep $line $in_data_drep | cut -d',' -f1,2,3,5,6 >> $summary_table_tmp
done < <(cut -d',' -f1 $in_data_checkm)
sort -t',' -k2,2nr -k3,3n -k5,5nr $summary_table_tmp > $mag_table
rm $summary_table_tmp
#All MAGs
echo '
MAG SUMMARY
        Total # MAGs' > $summary_table
cat $mag_table | wc -l >> $summary_table
#High quality
echo '       High quality' >> $summary_table
awk -F ',' '($2 > 90) && ($3 < 5) { print}' $mag_table | wc -l >> $summary_table
#Good quality
echo '  Good quality' >> $summary_table
awk -F ',' '($2 > 80) && ($3 < 10) { print}' $mag_table | wc -l >> $summary_table
#Near complete
echo '  Nearly complete' >> $summary_table
awk -F ',' '($2 > 98) && ($3 < 5) { print}' $mag_table | wc -l >> $summary_table
