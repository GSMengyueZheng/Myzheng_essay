Dat <- read.table("TMB_upper_25.txt",header = T,stringsAsFactors = F)
library(ggplot2)
library(reshape2)
Datme <- melt(Dat)
Datme$Signatures <- factor(Datme$SampleName, levels = paste0("Signature.",c(1:30)))

Datme$Signatures2[Datme$Signatures == "Signature.1"] <- "aging"
Datme$Signatures2[Datme$Signatures == "Signature.2"] <- "APOBEC"
Datme$Signatures2[Datme$Signatures == "Signature.3"] <- "BRCA"
Datme$Signatures2[Datme$Signatures == "Signature.4"] <- "smoking"
Datme$Signatures2[Datme$Signatures == "Signature.5"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.6"] <- "dMMR"
Datme$Signatures2[Datme$Signatures == "Signature.7"] <- "UV"
Datme$Signatures2[Datme$Signatures == "Signature.8"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.9"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.10"] <- "POLE"
Datme$Signatures2[Datme$Signatures == "Signature.11"] <- "alkylating"
Datme$Signatures2[Datme$Signatures == "Signature.12"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.13"] <- "APOBEC"
Datme$Signatures2[Datme$Signatures == "Signature.14"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.15"] <- "dMMR"
Datme$Signatures2[Datme$Signatures == "Signature.16"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.17"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.18"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.19"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.20"] <- "dMMR"
Datme$Signatures2[Datme$Signatures == "Signature.21"] <- "dMMR"
Datme$Signatures2[Datme$Signatures == "Signature.22"] <- "nucleotide excision repair"
Datme$Signatures2[Datme$Signatures == "Signature.23"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.24"] <- "aflatoxin"
Datme$Signatures2[Datme$Signatures == "Signature.25"] <- "Hodgkin lymphoma"
Datme$Signatures2[Datme$Signatures == "Signature.26"] <- "dMMR"
Datme$Signatures2[Datme$Signatures == "Signature.27"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.28"] <- "other"
Datme$Signatures2[Datme$Signatures == "Signature.29"] <- "tobacco chewing"
Datme$Signatures2[Datme$Signatures == "Signature.30"] <- "other"


pdf("TMB_upper_25.txt.pdf",10,4)
p <- ggplot(Datme,aes(x=variable,y=value))+xlab("Sample")+ylab("Proportion")+geom_bar(aes(fill= Signatures2),stat="identity", position="stack")
p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+scale_fill_discrete(name="Signatures")
dev.off()

