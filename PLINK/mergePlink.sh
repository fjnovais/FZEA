GENO1=$1     #plink
GENO2=$2     #plink 
GENO3=$3
BREEDHD1=$4
BREEDHD2=$5
BREEDHD3=$6

echo "${GENO2}.bed ${GENO2}.bim ${GENO2}.fam
${GENO3}.bed ${GENO3}.bim ${GENO3}.fam" > allfiles.txt 

./plink --cow --noweb --nonfounders --no-sex --bfile $GENO1 --merge-list allfiles.txt --make-bed --out mergeall

./plink --cow --noweb --nonfounders --bfile mergeall --no-parents --no-fid --no-sex --no-pheno --maf 0.05 --mind 0.1 --make-bed --out ${BREEDHD1}${BREEDHD2}${BREEDHD3}

./converttogeno.sh ${BREEDHD1}${BREEDHD2}${BREEDHD3} ${BREEDHD1}${BREEDHD2}${BREEDHD3}


rm mergeall.* allfiles.txt
