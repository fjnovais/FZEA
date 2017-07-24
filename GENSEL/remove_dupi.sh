awk '{if (x[$3]) { x_count[$3]++; print $0; if (x_count[$3] == 1) { print x[$3] } } x[$3] = $0}' snp_map | awk '{if($2<30 && substr($4,1,8)=="Internal") print $1}' > duplicates                       # position of the duplicates SNP that are named Internal

awk '{print $2}' 5_genoheifers.txt | sed -e 's/0/0 /g' -e 's/1/1 /g' -e 's/2/2 /g' -e 's/5/5 /g' > open_geno

awk '{$69562=$67707=$60564=$57132=$54145=$53172=$27558=$21027=$20916=$18117=$16165=$8977=$8900=$7479=$4216=$4171=$4117=$3966=$3941=$3885=$3680=""; print}' open_geno | sed 's/ //g' > tmp3

awk '{print $1}' 5_genoheifers.txt> tmp4

paste -d" " tmp4 tmp3 > 5_genoheifers_nodup.txt         #genotype without duplicated snps

awk 'FNR==NR{a[$1]++;next}!a[$1]' duplicates snp_map | awk '{print NR,$2,$3,$4}' >  snp_map_nodup                  #remove duplicates of snp_map

