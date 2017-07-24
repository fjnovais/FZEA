######LEMBRAR DE SEMPRE VERIFICAR OS ENDERECOS DOS ARQUIBOS###################

awk -F" " '{print $2}' ../1_Gift/2_QC_Preg/11_SECFMSM/41_geno0125_SECFMSM.txt_clean | sed -e 's/1/1 /g' -e 's/2/2 /g' -e 's/0/0 /g' -e 's/5/5 /g' > geno4.txt

awk -F" " '{print $1}' ../1_Gift/2_QC_Preg/11_SECFMSM/41_geno0125_SECFMSM.txt_clean  > id.txt

awk -F" " '{print $4}' ../1_Gift/2_QC_Preg/11_SECFMSM/3_snp_map_1558.txt_clean > snps

awk '{print $4,$2,$3}' ../1_Gift/2_QC_Preg/11_SECFMSM/3_snp_map_1558.txt_clean > tmp

echo "SNP_name SNP_chr SNP_pos" > header_map
cat header_map tmp > snp_map29_SECFMSM.txt                              #GENSEL MAP FILE
 
echo "SNP_name" > tmp2
cat tmp2 snps > snps2
awk '{for (f = 1; f <= NF; f++) a[NR, f] = $f} NF > nf { nf = NF } END {for (f = 1; f <= nf; f++) for (r = 1; r <= NR; r++) printf a[r, f] (r==NR ? RS : FS)}'  snps2 > snpst

############################################################
R < r_part.R --no-save
#########################################################################

cat snpst gensel.txt > gensel29_SECFMSE.txt                           #GENSEL GENOTIPO

rm  header_map tmp tmp2 snps snps2 snpst gensel.txt id.txt geno4.txt
