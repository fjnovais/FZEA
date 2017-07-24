hd_geno=$1
ld_geno=$2
out=$3

cat $ld_geno | awk '{print $1}' > id_ld
cat $ld_geno | awk '{$1=""; print}' | sed 's/ //g' | paste -d" " id_ld - | awk '{print $1,2,$2}' > ld

cat $hd_geno | awk '{print $1}' > id_hd  
cat $hd_geno | awk '{$1=""; print}' | sed 's/ //g' | paste -d" " id_hd - | awk '{print $1,1,$2}' > hd

echo "done :D"           

echo "ID Chip Call..." > header

cat header hd ld > $out.txt

rm header hd ld id_*
