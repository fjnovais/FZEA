#!/bin/bash

Truegeno=$1
outtruegeno=$2


###################################################################
# # Download PLINK from the web if not available
if [ ! -f plink ]; then
 echo "plink was not found in the current directory, thus it been download"
 echo " "
 wget http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-x86_64.zip
 unzip plink-1.07-x86_64.zip
 cp plink-1.07-x86_64/plink .
 ./plink --noweb --silent --file plink-1.07-x86_64/test
 rm -r plink-1.07-x86_64*
 if [ ! -f plink.log ]; then
  wget http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-i686.zip
  unzip plink-1.07-i686.zip
  cp plink-1.07-i686/plink .
  rm -r plink-1.07-i686*
  echo " "
 fi
rm plink.log
echo " "
fi

#### Checking if input file is available
if [ ! -f ${Truegeno}.bed ] || [ ! -f ${Truegeno}.bim ] || [ ! -f ${Truegeno}.fam ]; then
 echo "*** File " ${Truegeno}.* " representing the genotype file of PLINK binary format were not found ***"
 echo " "
 echo " "
 exit
fi


echo "****    Data processing for imputation started     ****"
Allelecode=$(awk '{print $6}' ${Truegeno}.bim | sort | uniq | awk '{if ($1==1) print "12"; else if ($1=="B") print "AB"; else if($1=="G" || $1=="T" || $1=="C") print "ACGT"}')

#####  REFERENCE  #######
if [ $Allelecode = 12 ]; then
 cat ${Truegeno}.bim | awk '{print $2,2}' > recodeallele.txt
 ./plink --silent --cow --noweb --nonfounders --bfile ${Truegeno} --make-bed --out ref_upd
 ./plink --silent --cow --noweb --nonfounders --bfile ref_upd --recodeA --recode-allele recodeallele.txt --out geno
rm recodeallele.txt

elif [ $Allelecode = AB ]; then
 cat ${Truegeno}.bim | awk '{print $2,"A","B",1,2}' > alleleupdate.txt
 ./plink --silent --cow --noweb --nonfounders --bfile ${Truegeno} --update-alleles alleleupdate.txt --make-bed --out ref_upd
 cat ref_upd.bim | awk '{print $2,2}' > recodeallele.txt
 ./plink --silent --cow --noweb --nonfounders --bfile ref_upd --recodeA --recode-allele recodeallele.txt --out geno
rm alleleupdate.txt recodeallele.txt

elif [ $Allelecode = ACGT ]; then
 echo "ACGT format is not allowed -- Please use an AB coding or 12"
 echo "The 12 format can be obtained with PLINK using the --recode12"
 echo " "
 echo " "
 exit
fi

awk 'NR>1 {print $2}' geno.raw > IDs_sons.txt
awk 'NR>1' geno.raw | cut -d' ' -f7- | awk '{gsub(/NA/,5); print}' |
paste -d' ' IDs_sons.txt - > ${outtruegeno}.geno

rm geno.* IDs_sons.txt ref_upd.*

