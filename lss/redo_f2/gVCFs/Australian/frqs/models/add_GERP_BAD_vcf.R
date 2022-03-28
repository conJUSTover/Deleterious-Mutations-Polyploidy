library(stringr)
library(dplyr)
library(reshape2)

setwd("/work/LAS/jfw-lab/del_mut_gvcf/redo_f2/gVCFs/Australian/frqs/models")

infile <- read.table("pSONIC.NO_GC.recode.noheader.vcf", header = TRUE)
names(infile) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","100","101","102","103","104","105","106","107","108","109","301","302","303","304","305","306","307","401","402","403","405","406","502","504","505","507","601","602","603","604","605","606","701","702","800","801","802","803","900","901","902","903","904","C2-4","G2-21","G3-83")


infile$CHROM <- str_pad(infile$CHROM, 2, pad = "0")
infile$POS <- as.character(infile$POS)
infile$CHROM_POS <-  apply(infile[,c("CHROM", "POS")], 1, paste, collapse = "_")
#POS <- read.table("POS.GERP.Pvalue.txt", header = TRUE)
infile <- subset(infile, select = -c(ID, REF, ALT, QUAL, FILTER, INFO, FORMAT))
infile$POS <- as.integer(infile$POS)



RefEqual <- infile[c(infile$`C2-4` == "0/0" | infile$`G2-21` == "0/0"),]
AltEqual <- infile[c(infile$`C2-4` == "1/1" | infile$`G2-21` == "1/1"),]

#sub 0/0 for 2 if alternative == ancestral
AltEqual[3:48] <- lapply(AltEqual[,3:48], factor, levels=c('./.', '0/0', '0/1', '1/1'), labels=c(NA, 2, 1, 0))


#sub 1/1 for  2 if reference == ancestral
RefEqual[3:48] <- lapply(RefEqual[,3:48], factor, levels=c('./.', '0/0', '0/1', '1/1'), labels=c(NA, 0, 1, 2))

#recombine dfs
infile <- rbind(AltEqual, RefEqual)

#Split into subgenomes, filter for NAs
#Dfile <- infile[infile$CHROM > 13,]
#Dfile <- Dfile[complete.cases(Dfile[,c(3:40)]),] #37, 41

#Afile <- infile[infile$CHROM < 14,]
#Afile <- Afile[complete.cases(Afile[,c(3:36,41:45)]),]

#infile <- rbind(Afile, Dfile)
#rm(Afile)
#rm(Dfile)

GERP_vcf <- infile[FALSE,]


for(i in c(01:26)){
  j <- str_pad(i, 2, pad = "0")
  print(j)
  GERP_name <- paste("../../../GERP_SCORES/roast", j, "msa.in.new.raw.rates", sep = ".")
  GERP <- read.table(GERP_name, header=FALSE)
  names(GERP) <- c('Neutral Rate','RS Score')
  temp_frqs <- infile[infile$CHROM == j,]
  temp_frqs$GERP <- GERP[temp_frqs$POS,2]
  GERP_vcf <- dplyr::bind_rows(temp_frqs, GERP_vcf)
}

print("GERP done")
rm(infile)

BAD_results <- read.table("../all.GT.merged_results.txt", header = TRUE, sep = "\t")
BAD_results <- as.data.frame(BAD_results[,c(21, 15, 14, 3)])


GERP_vcf <- as.data.frame(GERP_vcf)

almost_vcf <- merge(x=GERP_vcf, y=BAD_results, by.x = c("CHROM_POS"), by.y = c("SNP_Name"), all.x = TRUE)

print("merged_results done")

rm(GERP_vcf)
rm(BAD_results)

synonymous <- read.table("../all.GT.se_variant.txt", header = FALSE)
names(synonymous) <- c("CHROM_POS", "Synonymous", "Syn_Gene")
final_vcf <- merge(x=almost_vcf, y=synonymous, by.x = c("CHROM_POS"), by.y = c("CHROM_POS"), all.x = TRUE)

print("se_variant done")

write.table(final_vcf, file = "pSONIC.GT.NO_GC.BAD.allGERP.vcf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

GERP06 <- subset(final_vcf, GERP > 0 | !is.na(GeneID) | !is.na(Synonymous))
write.table(GERP06, file = "pSONIC.GT.NO_GC.GERP06_BAD.vcf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

