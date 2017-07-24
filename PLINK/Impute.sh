GENO=$1
MAP=$2
PED=$3
LD_DEN=$4
HD_DEN=$5
BREEDHD=$6
BREEDLD=$7

#####################

echo "title=*Imputation $LD_DEN to $HD_DEN*;
genotype_file=*$GENO*;
snp_info_file=*$MAP*;
ped_file=*$PED*;
parentage_test /ert_mm=0.001 /ert_m=0.0008 /find_match_cnflt /find_match_mp /find_match_ugp /remove_conflict;
output_folder=*${LD_DEN}to${HD_DEN}_${BREEDLD}${BREEDHD}*;
save_genotype;
save_hap_lib /diplotype;
njob=4; " > card_${LD_DEN}to${HD_DEN}_${BREEDHD}${BREEDLD}.ctr

sed -i 's/*/"/g' card_${LD_DEN}to${HD_DEN}_${BREEDHD}${BREEDLD}.ctr

FImpute card_${LD_DEN}to${HD_DEN}_${BREEDHD}${BREEDLD}.ctr -o
