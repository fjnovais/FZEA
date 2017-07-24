geno=$1

cat $geno | awk '{print $2}' | sed -e 's/2/7/g' -e 's/0/2/g' -e 's/7/0/g' > tmp1

cat $geno | awk '{print $1}' > tmp2 

paste -d" " tmp2 tmp1 > ${geno}_inverted

echo ":)"

rm tmp*
