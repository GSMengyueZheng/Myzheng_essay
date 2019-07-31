##################################
####################²»Çø·Ö°©Ö×
Dat <- read.table("All_noindel.txt.sta.txt",stringsAsFactors = F,header=F)
library(ggplot2)
library(stringr)
f <- function(a,b){
  if(a == "C" && b == "A"){
    "C>A"
  }else if(a == "G" && b == "T"){
    "C>A"
  }else if(a == "C" && b == "G"){
    "C>G"
  }else if(a == "G" && b == "C"){
    "C>G"
  }else if(a == "C" && b == "T"){
    "C>T"    
  }else if(a == "G" && b == "A"){
    "C>T"
  }else if(a == "T" && b == "A"){
    "T>A"
  }else if(a == "A" && b == "T"){
    "T>A"
  }else if(a == "T" && b == "C"){
    "T>C"
  }else if(a == "A" && b == "G"){
    "T>C"
  }else if(a == "T" && b == "G"){
    "T>G"
  }else{
    "T>G"
  }
}

l <- length((Dat$V2))
for (i in 1:l){
  Be <- str_sub(Dat$V1[i],2,2)
  Af <- str_sub(Dat$V1[i],7,7)
  Dat$V3[i] <- f(Be,Af)
  Dat$V4[i] <- str_sub(Dat$V1[i],1,3)
}
Dat2 <- Dat[order(Dat$V3),]
CRCsum <- sum(Dat2$V2)
for (i in 1:length(Dat2$V1)){
  Dat2$V6[i] = round((Dat2$V2[i]/CRCsum*100),2)
 }
Dat2$V1 <- factor(Dat2$V1,levels=Dat2$V1)
pdf("All_Cancer.pdf",12,4)
ggplot(Dat2,aes(x=V1,y=V6,fill=V3))+geom_bar(stat="identity")+ylim(c(0,8))+xlab("")+
  ylab("Percent of Mutations(%)")+scale_x_discrete(labels=Dat2$V4)+
  theme(axis.text.x = element_text(size = 8,angle = 90,face="bold"))+
  scale_fill_discrete("Substitution\n      Type")+theme(text=element_text(size=12))
dev.off()