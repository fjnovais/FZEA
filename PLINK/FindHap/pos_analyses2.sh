#!/bin/bash

chr=$1              #chromosome  
start=$2            #start(bp) of the interested region
end=$3              #end(bp) of the interested region
ped=$3              #pedigree file

echo "Your region was ${chr}:${start}:${end}, and the pedigree file ${ped}"

#begging of the region. $4=overall number. Look for this number in hap.list file $7(overall)
sed 's/ \+/ /g' chromosome.data | awk '{if($2=='$chr' && $5=='$start') print $4}' > start_overall
echo "The firt overall number is" | cat - start_overall
sed 's/ \+/ /g' chromosome.data | awk '{if($2=='$chr' && $5=='$end') print $4}' > end_overall
echo "The last overall number is" | cat - end_overall

start_over=`cat start_overall`
end_over=`cat end_overall`
declare -i start_over2
declare -i end_over2
start_over2=$start_over-500
end_over2=$end_over+500

#position from the map_reg. Look column $7(seg)
echo "Looking values of column 7 (seg)"
echo "seg > = $start_over2 e seg < = $end_over2";
sed -e 's/^ *//g' -e 's/ \+/ /g' hap.list | awk '{if($4=="'$chr'" && $7>="'$start_over2'" && $7<"'$end_over2'") print $4,$7,$6}' | sort -u > tmp
echo "chr locus seg" | cat - tmp 

read -p "Which segment is your haplotype?: " hap

declare -i hap_column
hap_column=$hap+4
#the + 4 it is because we have 4 columns before the haplotypes.

awk 'NR>4 {print $1,$2,$'$hap_column'}' hap.found > hap_tmp
awk 'NR>9 {print $1,$4}' $ped > ped_tmp

awk 'NR==FNR{a[$1]=$2;next} {print $0,a[$1]}' ped_tmp hap_tmp | awk '{print $4,$1,$2,$3}' > hap_ped

echo "id id_r pnt hap"$hap > head_hap

cat head_hap hap_ped > hap_chr${chr}_pos${start}_${end}

rm head_hap hap_tmp hap_ped start_overall end_overall tmp ped_tmp

echo "Finished!"
