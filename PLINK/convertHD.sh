GENO1=$1     #.csv
MAP=$2      
BREEDHD=$3

cat $GENO1 | awk -F"," 'NR>1 {print $1}' > ids.txt
cat $GENO1 | awk -F"," 'NR>1 {$1=""; print}' | sed 's/?_?/0_0/g' | sed 's/_/ /g' | paste -d" " ids.txt - | awk '!_[$1]++' > file.ped

cat $MAP | awk -F"," 'NR>1 {print $2,$1,0,$3}' > file.map
cat file.map | awk '{if($1=="0" || $1=="X" || $1=="Y" || $1=="MT") print}' > snplist
cat file.map | awk '{if($1!="0" && $1!="X" && $1!="Y" && $1!="MT") print $1"_"$4,$2}' | awk '_[$1]++ {print $2}' >> snplist

./plink --cow --noweb --nonfounders --file file --no-parents --no-fid --no-sex --no-pheno --exclude snplist --make-bed --out ${BREEDHD}_raw

./plink --cow --noweb --nonfounders --bfile ${BREEDHD}_raw --no-parents --no-fid --no-sex --no-pheno --maf 0.01 --mind 0.1 --hwe 0.00001 --make-bed --out $BREEDHD
./converttogeno.sh $BREEDHD $BREEDHD

./plink --cow --noweb --nonfounders --bfile ${BREEDHD}_raw --no-parents --no-fid --no-sex --no-pheno --keep 300Holstein.txt --make-bed 300${BREEDHD}_raw


rm snplist ids.txt file.*


