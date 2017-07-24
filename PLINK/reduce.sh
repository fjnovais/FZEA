##Impute LD Map##
MAP=$1
BREED=$2
DENSITY=$3

cat $MAP | awk '{print $1}' > snplist

./plink --cow --noweb --nonfounders --bfile $BREED --no-parents --no-fid --no-sex --no-pheno --extract snplist  --make-bed --out ${DENSITY}_${BREED}

./converttogeno.sh ${DENSITY}_${BREED} ${DENSITY}_${BREED}

rm snplist 
