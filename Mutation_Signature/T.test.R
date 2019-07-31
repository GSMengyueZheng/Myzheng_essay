library(ggplot2)
library(reshape2)

Dat <- read.csv("breast.list.res",header = F,sep="\t")
Dat$aflatoxin <- Dat$V27
Dat$aging <- Dat$V4
Dat$alkylating <- Dat$V14
Dat$APOBEC <- Dat$V5+Dat$V16
Dat$BRCA <- Dat$V6
Dat$dMMR <- Dat$V9+Dat$V18+Dat$V23+Dat$V24+Dat$V29
Dat$Hodgkin_lymphoma <- Dat$V28
Dat$nucleotide_excision_repair <- Dat$V25
Dat$other <- Dat$V8+Dat$V11+Dat$V12+Dat$V15+Dat$V17+Dat$V19+Dat$V20+Dat$V21+
  Dat$V22+Dat$V26+Dat$V30+Dat$V31+Dat$V33
Dat$POLE <- Dat$V13
Dat$smoking <- Dat$V7
Dat$tobacco_chewing <- Dat$V32
Dat$UV <- Dat$V10
DatP <- Dat[Dat$V1 == "Position",c(1:3,34:46)]
DatN <- Dat[Dat$V1 == "Negative",c(1:3,34:46)]

lst <- c("aflatoxin","aging","alkylating","APOBEC","BRCA","dMMR","Hodgkin_lymphoma",
         "nucleotide_excision_repair","other","POLE","smoking","tobacco_chewing","UV")
Mat <- NULL

for (i in 1:length(lst)){
MeanP <- mean(DatP[,i+3])
MeanN <- mean(DatN[,i+3])
Res <- t.test(DatP[,i+3],DatN[,i+3],alternative = "two.sided")
ResP <- Res$p.value
Line <- c(lst[i],MeanP,MeanN,ResP)
Mat <<- rbind(Mat,Line)
}

write.table(Mat,"breast.list.ttest.csv",sep=",",col.names=c("Signature","Mean_Positive","Mean_Negative","P.value"),row.names=F)
