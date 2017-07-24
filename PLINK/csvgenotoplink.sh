GENO=$1
MAP=$2
BREED=$3

cat $GENO | awk -F"," 'NR>1 {print $1}' > ${BREED}_ids.txt
cat $GENO | awk -F"," 'NR>1 {$1=""; print}' | sed 's/?_?/0_0/g' | sed 's/_/ /g' | paste -d" "  ${BREED}_ids.txt - | awk '!_[$1]++' > ${BREED}.ped

cat $MAP | awk -F"," 'NR>1 {print $2,$1,0,$3}' > ${BREED}.map

cat $BREED.map | awk '{if($1=="0" || $1=="X" || $1=="Y" || $1=="MT") print}' > snplist

./plink --cow --noweb --nonfounders --file $BREED --no-parents --no-fid --no-sex --no-pheno --exclude snplist --make-bed --out $BREED
./plink --cow --noweb --nonfounders --bfile $BREED --no-parents --no-fid --no-sex --no-pheno --recodeA --out $BREED

cat $BREED.raw | awk 'NR>1 {print $1}' > tmp1

cat $BREED.raw | awk 'NR>1 {$1=$2=$3=$4=$5=$6=""; print}' | sed 's/NA/5/g' | sed 's/ //g' | paste -d" " tmp1 - | awk '{print $1,1,$2}' > ${BREED}_geno.txt

rm tmp1 snplist
