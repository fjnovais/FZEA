echo "##################Preparing the FImpute genotype file#####"

geno_sires=$1
geno_prog=$2
ped_r=$3
map=$4


cat $geno_sires $geno_prog | awk '{print $1,1,length($3),$3}' | sort -k1 > tmp

cat $ped_r | awk 'NR>9 {print $4,$1}' | sort -k1 > tmp2

#relaciona os dois bancos, pegando so as 3 primeiras colunas (id, animal, n_SNPs)
awk 'NR==FNR{a[$1]=$3;next} ($1) in a {print $1,$2}' tmp tmp2 > tmp3
sort tmp3 > tmp4
join tmp4 tmp > tmp5

awk '{printf ("%10s %9s %9s %s\n"),$2,$3,$4,$5}' tmp5 > genotypes.txt

rm tmp*

echo "##################Making Pedigree File##################"
echo "##################Pedigree must be renum###############"

cat $geno_prog | awk '{print $1}' | sort > tmp            
cat $ped_r | sed '1,9d' | sort -k4 > tmp2       
join -1 1 -2 4 tmp tmp2 > tmp3

awk '{if($3>0) print "M",$3,"-2","-1","20050101",$3,$3}' tmp3 | sort -k2 -u > tmp4
awk '{if($4>0) print "F",$4,"-1","-2","20050101",$4,$4}' tmp3 | sort -k2 -u > tmp5
awk '{$1=""; print "F",$0,"20120101",$2,$2}' tmp3 > tmp6

cat tmp4 tmp5 tmp6 > tmp7
awk '{printf "%s %11s %8s %8s %8s %16s %29s\n",$1,$2,$3,$4,$5,$6,$7}' tmp7 > pedigree.file

rm tmp*

echo "##################Making map file##################"

#loop for add number column for each chr
NM=$(awk 'NR==2{print length($3)}' $geno_sires)

cat $map | awk 'NR>1 {print $4,$2,$3 }' > tmp0

for i in {1..29}
do
grep " $i " tmp0 | cat -n >> map_tmp
done

#loop para formar a coluna de Markers
for i in {1..9}
do
echo "Marker00000000"$i >> marker
done
for i in {10..99}
do
echo "Marker0000000"$i >> marker
done
for i in {100..999}
do
echo "Marker000000"$i >> marker
done
for i in {1000..9999}
do
echo "Marker00000"$i >> marker
done
for ((i=10000; i<=$NM; i++)); do echo "Marker0000"$i >> marker; done

paste marker map_tmp > map_tmp2

awk '{printf ("%s %4s %8s %8s %11s %4s %8s\n"), $1,$4,$2,NR,$5,"1",NR}' map_tmp2 > map_tmp3
echo "SNPname          chrome  within  overall  location  n_chips   chip1" > tmp
cat tmp map_tmp3 > chromosome.data

rm tmp map_tmp map_tmp2 map_tmp3 marker
