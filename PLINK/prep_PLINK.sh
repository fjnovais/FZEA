genohd=$1
genold=$2
maphd=$3
mapld=$4
outgenoREF=$5
outgenoVAL=$6


cat $genohd | awk '{print $1}' > ids.txt

cat $genohd | awk '{print $2}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print }' | 
./change012to111222.awk | paste -d" " ids.txt - > ${outgenoREF}.ped

cat $maphd | awk '{print $2,$1,0,$3}' > ${outgenoREF}.map 

rm ids.txt

#### validation animals
cat $genold | awk '{print $1}' > ids.txt

cat $genold | awk '{print $2}' | 
awk 'BEGIN {FS="";OFS=" "} {$1=$1; print }' | 
./change012to111222.awk | paste -d" " ids.txt - > ${outgenoVAL}.ped

cat $mapld | awk '{print $2,$1,0,$3}' > ${outgenoVAL}.map 

rm ids.txt
