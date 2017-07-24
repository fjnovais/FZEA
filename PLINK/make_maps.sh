###Input Map files##
HD=$1                     
LD=$2
DENLD=$3
DENHD=$4
FOLDER=$(pwd)

echo "Name Chr Pos" > head_map
 
cat $HD | awk '{print $2,$1,$4}' | cat head_map - > ${DENHD}_map.txt
cat $LD | awk '{print $2,$1,$4}' | cat  head_map - > ${DENLD}_map.txt

rm head_map

####################MAKING MAPS#################

echo "// #################################################################################
// ================= MANDATORY INFORMATION =========================================
GENOFILE_0125 heifers_imputed.txt
NPARALELJOBS 3
MAP(NCP) snp_info.txt
FOLDERLOCATION $FOLDER
FLAGWINLINUX /
// ============== QUALITY CONTROL (NO CR_SNP YET) ==================================
CR_ANI NO
MAF 0.05
CR_SNP 0.95
// =================== GENERAL FUNCTIONS ===========================================
FLAGPARENTDISCOV YES
SELECANI NO
SNPINFOFLAGREDUCINGFORMAT snp_info_54609_flag6K.txt
REDUCEGENOFILE genotype_9.txt
FILEANIMALSCHIPFIMPUTE list_REFIMP.txt
NCLUSTERS_K 3
SNPLISTTOPRINT snp_list.txt
// =================== IMPUTATION ACCURACY =========================================
FOLDERLOCATIONORIGFILE /home/gerson_jr/2_embrapa/1_genotipos/1_Gir
GENO(FIMPUTE)FILE BULLTRUE_nospace.geno18
FOLDERLOCATIONIMPUTEDFILE /home/gerson_jr/2_embrapa/1_genotipos/1_Gir
GENO(IMPUTED)FILE BULLUNGENO_ped0.01_nospace.geno18
MAPFIMPUTEHD_BEFORE tmp1
MAPFIMPUTEHD_AFTER tmp1
ANIMALLIST listAnimals_imputed.txt
// ================== CONVERTING LAB FILE ==========================================
LABRAWFILE 25Jul14_Bio015-Yuri_50k_SR.txt
SNPMAPRAW snp_info_54609.txt
NCHRRAWDATA 31
LINESTOJUMP 10
COLUMNOPTION 7
// ================== ROH ==========================================================
WINSIZE 50
DISTBETWEENSNPS NO
NHET 1
NMISSING 2
// ================== COMBINING AFFY ILLUNIMA ======================================
SNPINFOAFFY snp_infoAffy.txt
GENOAFFY genoAffy.txt
SNPINFOILLU snp_infoIllu.txt
GENOILLU genoIllu.txt
// ================== COMMERCIAL IMPUTATION ========================================
SNPINFOLD ${DENLD}_map.txt
GENOLD 70k_Gir_geno.txt
SNPINFOHD ${DENHD}_map.txt
GENOHD HD_Gir_geno_nodup.txt
ANIMALLIST_RI list_RI.txt
// (PLINK/REDUCE/REDUCE_U/FIMPUTECHIP/PLOTR/FINDSNP/IMPACC/EXPAND/SORT/MAF/MAFFI) ==
// (HAM_LONG) ======================================================================
SPECIFICFUNCTION FIMPMAPFROM2SNPINFO
// #################################################################################" > param.txt

./RTBM_1001_V056

mv snp_info.txt ${DENHD}_${DENLD}.map

rm LOG.txt
