// #################################################################################
// ================= MANDATORY INFORMATION =========================================
GENOFILE_0125 heifers_imputed.txt
NPARALELJOBS 5
MAP(NCP) snp_info.txt
FOLDERLOCATION /home/1_Gift/6_Imputing/Prog_Ricardo
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
FOLDERLOCATIONORIGFILE /home/1_Gift/6_Imputing/Prog_Ricardo
GENO(FIMPUTE)FILE BULLTRUE_nospace.geno18
FOLDERLOCATIONIMPUTEDFILE /home/1_Gift/6_Imputing/Prog_Ricardo
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
SNPINFOLD snp_map_ld.txt
GENOLD ~/1_Gift/2_QC_Preg/8_ENSECFMSM/3_all/5_genoheifers_nodup.txt_clean
SNPINFOHD snp_map_hd.txt
GENOHD ~/1_Gift/2_QC_Preg/12_sires/geno_sires.txt_clean
ANIMALLIST_RI list_RI.txt
// (PLINK/REDUCE/REDUCE_U/FIMPUTECHIP/PLOTR/FINDSNP/IMPACC/EXPAND/SORT/MAF/MAFFI) ==
// (HAM_LONG) ======================================================================
SPECIFICFUNCTION PLOTACCFIMPUTE_SNP_ANI
// #################################################################################

