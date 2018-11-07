perl neoantigen_v2.pl HLA.txt input.txt neo8.txt 1 2 3 10 11 12 8 > neo8.log 2>&1
perl neoantigen_v2.pl  HLA.txt input.txt neo9.txt 1 2 3 10 11 12 9 > neo9.log 2>&1
head -1 neo8.txt > head
sed 1d neo9.txt |cat neo8.txt - > neo_all.txt
perl getNE.pl neo_all.txt |sort -k 1,1 -k 2,2 -k5,5 -k 8,8 > neo_tmp
cat head neo_tmp > neo_filter.txt

#perl neoantigen2_v2.pl HLA2.txt input.txt neo_MHCII.txt 1 2 3 10 11 12 14 > neo_mhc2.log 2>&1
#perl getNE.pl neo_MHCII.txt |sort -k 1,1 -k 2,2 -k 5,5 -k 8,8 > neo_MHCII.tmp
#cat head neo_MHCII.tmp > neo_MHCII_filter.txt

#rm neo_tmp neo_MHCII.tmp
