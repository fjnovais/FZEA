##Impute LD Map##
MAP=$1
BREED=$2
DENSITY=$3

cat $MAP | awk '{print $1}' > snplist

./plink --cow --noweb --nonfounders --bfile $BREED --no-parents --no-fid --no-sex --no-pheno --extract snplist  --make-bed --out ${DENSITY}_${BREED}

./converttogeno.sh ${DENSITY}_${BREED} ${DENSITY}_${BREED}

rm snplist 

__________________________________________________________________________________________________________________________________________________________________________
at $MAP | awk -F"," 'NR>1 {if($2!="0" && $2!="X" && $2!="Y" && $2!="MT") print $1}' | sed -e 's/PAGLDMax_//g' -e 's/ZTSLD2.0_//g' > snplist_${DENSITY}

./plink --cow --noweb --nonfounders --bfile $BREED --no-parents --no-fid --no-sex --no-pheno --extract snplist_$DENSITY --make-bed --out $DENSITY'_'$BREED
./plink --cow --noweb --nonfounders --bfile $DENSITY'_'$BREED --no-parents --no-fid --no-sex --no-pheno --recodeA --out $DENSITY'_'$BREED

cat $DENSITY'_'$BREED.raw | awk 'NR>1 {print $1}' > tmp1

cat $DENSITY'_'$BREED.raw | awk 'NR>1 {$1=$2=$3=$4=$5=$6=""; print}' | sed -e 's/NA/5/g' -e 's/ //g' | paste -d" " tmp1 - | awk '{print $1,2,$2}' > ${DENSITY}_${BREED}_geno.txt

rm tmp1
