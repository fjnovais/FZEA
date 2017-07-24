FOLDER=$1
DENLD=$2        #20k,50k..
DENHD=$3
TRUE=$4


##TRUE GENO (CLEAN)##
awk 'NR==FNR {a[$1]=$2;next} !($2 in a) {print $2}' ${FOLDER}/snp_info.txt ${TRUE}.bim > snplist

./plink --cow --noweb --nonfounders --bfile $TRUE --no-parents --no-fid --no-sex --no-pheno --exclude snplist --make-bed --out true_geno${DENHD}

./converttogeno.sh true_geno${DENHD} true_geno${DENHD}

rm snplist

##############IMPUTED###############
cat ${FOLDER}/genotypes_imp.txt | awk 'NR>1 {if($2=="2") print $3}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print $0}' | 
awk '{for (i=1;i<=NF;i++) { if($i==0) $i="1 1"; else if($i==1) $i="1 2"; else if($i==2) $i="2 2"; else if($i==5) $i="0 0"}  print}' > geno

cat ${FOLDER}/genotypes_imp.txt | awk 'NR>1 {if($2=="2") print $1,$1,0,0,0,-9}' > ids
paste -d' ' ids geno > file.ped

cat ${FOLDER}/snp_info.txt | awk 'NR>1 {print $2,$1,0,$3}' > file.map
rm geno ids

awk 'NR==FNR {a[$2]=$2;next} ($1 in a) {print $1}' true_geno${DENHD}.bim ${FOLDER}/snp_info.txt > snplist2

./plink --cow --noweb --nonfounders --file file --extract snplist2 --make-bed --out geno${DENLD}_imp

./converttogeno.sh geno${DENLD}_imp geno${DENLD}_imp

mv geno${DENLD}_imp.* ${FOLDER}/
mv true_geno${DENHD}.* ${FOLDER}/

rm file.map file.ped snplist2

