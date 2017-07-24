ped_r=$1
chr=$2
OUT=$3

echo "Considerando o ped como id,sire,dam,id_r"

awk 'NR==FNR {a[$1]=$4;next} ($1) in a {print a[$1],$4}' $ped_r haplotypes.txt > haplotypes_id.txt

awk '{print $2}' haplotypes_id.txt | ./change12tospace12.awk > tmp
awk '{print $1}' haplotypes_id.txt | paste -d" " - tmp > ${OUT}.ped
awk 'NR>1 {print $2,$1,"0",$3}' snp_info.txt > ${OUT}.map

rm tmp haplotypes_id.txt

./plink --cow --noweb --nonfounders --file $OUT --no-parents --no-fid --no-sex --no-pheno --chr $chr --recodeHV
