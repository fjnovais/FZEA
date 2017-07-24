rm(list=ls())
setwd("/home/gerson_jr/GenSel_Convert")
geno <- read.table("geno4.txt", header=F)
id=read.table("id.txt", header=F)

#geno[1:10,1:10]

gensel=(geno-1)*10
change=function(gensel){
  gensel[gensel==40]=NA
  gensel[is.na(gensel)]=mean(gensel,na.rm=T)
  gensel
}
genonew=apply(gensel,2,change)          #aplica a funcao criada por colunas (2). Se fosse por linhas seria (1)
genonew2=round(genonew, 0)                   #remove decimals
genonew3=data.frame(genonew2)
gensel=cbind(id,genonew3)

#gensel[1:10,1:10]

rm(genonew, genonew2, genonew3)

write.table(gensel, "gensel.txt", sep=" ", col.names = FALSE, row.names = FALSE, quote=FALSE)

