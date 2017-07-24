##Input Files##
LD=$1              #PLINK LD
DENLD=$2           #density LD  (eg. 20k, 50k, 70k)  
BREEDHD=$3       

#Common SNPs##
awk 'NR==FNR {a[$2]=$2;next} $2 in a {print $2}' ${BREEDHD}.bim ${LD}.bim > snplist
N_LD=$(awk 'END {print NR}' ${LD}.fam)

./plink --cow --noweb --nonfounders --bfile $LD --no-parents --no-fid --no-sex --no-pheno --extract snplist --make-bed --out Girol${N_LD}_${DENLD}

./converttogeno.sh Girol${N_LD}_${DENLD} Girol${N_LD}_${DENLD}

##Make FImput##
awk '{print $1,1}' ${BREEDHD}.geno > id1
awk '{$1=""; print}' ${BREEDHD}.geno | sed 's/ //g' | paste -d' ' id1 - > tmp

awk '{print $1,2}' Girol${N_LD}_${DENLD}.geno > id2                       
awk '{$1=""; print}' Girol${N_LD}_${DENLD}.geno | sed 's/ //g' | paste -d' ' id2 - > tmpLD

#number of animals#
N_HD=$(awk 'END {print NR}' tmp)                          

#Header#
echo "ID Chip Call..." > head_geno

#Merge Files#
if [ ! -d 1_Imputations ]; then
mkdir 1_Imputations
fi

cat head_geno tmp tmpLD > 1_Imputations/${BREEDHD}${N_HD}'_'${DENLD}Girol${N_LD}.geno

##Remove Files##
rm head_geno tmp tmpLD snplist id1 id2
