REF=$1
VAL=$2

if [ ! -f plink ]; then
wget https://www.cog-genomics.org/static/bin/plink150805/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip
fi

./plink --cow --noweb --nonfounders --file $REF --no-parents --no-fid --no-sex --no-pheno --exclude snp_error.txt --make-bed --out $REF 

./plink --cow --noweb --nonfounders --file $VAL --no-parents --no-fid --no-sex --no-pheno --make-bed --out $VAL

./change_snps_name.sh ${REF}.bim
./plink --cow --noweb --nonfounders --bfile $REF --bmerge ${VAL}.bed ${VAL}.bim ${VAL}.fam --make-bed --out tmp

### keep ids of reference and validation animals
cat ${REF}.fam | awk '{print $1,$2}' > ${REF}.ids
cat ${VAL}.fam | awk '{print $1,$2}' > ${VAL}.ids
cat ${VAL}.bim | awk '{print $2}' > snps.txt

### extract them from merge file
./plink --cow --noweb --nonfounders --bfile tmp --keep ${REF}.ids --extract snps.txt --make-bed --out ${REF}_finalFImp
./plink --cow --noweb --nonfounders --bfile tmp --keep ${VAL}.ids --extract snps.txt --make-bed --out ${VAL}_finalFImp

rm tmp* snps.txt




