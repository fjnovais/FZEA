snp_map_ld=$1
snp_map_hd=$2

#Cria lista com SNPs com mesma posicao (duplicados)
rm snp_error snp_error.txt

awk '{print $1"_"$4,$2}' bulls.map > tmppos
awk '{print $1"_"$4}' bulls.map | uniq -d > posdups
awk 'NR==FNR{a[$1]=$2;next} ($1) in a {print $0, a[$1]}' tmppos posdups | awk '{print $2}' > snp_error

#Adiciona na lista Snps com mesmo nome mas posicao diferentes
awk 'NR==FNR{a[$1]=$3;next} ($1) in a {print $0, a[$1]}' $snp_map_ld $snp_map_hd > tmp
awk '{print $1,$3-$4}' tmp | awk '{if($2!=0) print $1}' >>  snp_error

#SNPs com mesma posicao mas nomes diferentes
awk '{print $1"_"$4,$2}' heifers.map > tmp_h
awk '{print $1"_"$4,$2}' bulls.map > tmp_b
awk 'NR==FNR{a[$1]=$2;next} ($1) in a {print $0, a[$1]}' tmp_h tmp_b > tmp_bh
awk '$2!=$3 {print $2}' tmp_bh >> snp_error

cat snp_error | sort -u > snp_error.txt


rm snp_error posdups tmp* 
