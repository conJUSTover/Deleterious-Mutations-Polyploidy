library(stringr)
library(dplyr)
library(reshape2)

setwd("/work/LAS/jfw-lab/del_mut_gvcf/BQ20_gVCFs/only_wild/wild_gvcf/frqs")
infile <- read.table(gzfile("../all.GT.vcf.gz"))
names(infile) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","100","101","102","103","104","105","106","107","108","109","301","302","303","304","305","306","307","401","402","403","405","406","502","504","505","507","601","602","603","604","605","606","701","702","800","801","802","803","900","901","902","903","904","C2-4","G2-21","G3-83")


infile$CHROM <- str_pad(infile$CHROM, 2, pad = "0")
#infile$POS <- as.character(infile$POS)
#infile$CHROM_POS <-  apply(infile[,c("CHROM", "POS")], 1, paste, collapse = "_")
#POS <- read.table("POS.GERP.Pvalue.txt", header = TRUE)
#infile <- subset(infile, select = -c(ID, REF, ALT, QUAL, FILTER, INFO, FORMAT))
#infile$POS <- as.integer(infile$POS)




RefEqual <- infile[c(infile$`C2-4` == "0/0" | infile$`G2-21` == "0/0"),]
AltEqual <- infile[c(infile$`C2-4` == "1/1" | infile$`G2-21` == "1/1"),]

#sub 0/0 for 2 if alternative == ancestral
AltEqual[10:55] <- lapply(AltEqual[,10:55], factor, levels=c('./.', '0/0', '0/1', '1/1'), labels=c('./.', '1/1', '0/1', '0/0'))


#don't change anything if reference == ancestral
#RefEqual[3:48] <- lapply(RefEqual[,3:48], factor, levels=c('./.', '0/0', '0/1', '1/1'), labels=c(NA, 0, 1, 2))

#recombine dfs
infile <- rbind(AltEqual, RefEqual)
infile <- infile[with(infile, order(CHROM, POS)),]

names(infile) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","100","101","102","103","104","105","106","107","108","109","301","302","303","304","305","306","307","401","402","403","405","406","502","504","505","507","601","602","603","604","605","606","701","702","800","801","802","803","900","901","902","903","904","C2-4","G2-21","G3-83")

write.table(infile, file="all.GT.polarized.vcf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
