library(stringr)
library(dplyr)
library(reshape2)

setwd("/work/LAS/jfw-lab/del_mut_gvcf/redo_f2/gVCFs/Australian/frqs/models")

in_vcf <- read.table("pSONIC.NO_GC.recode.noheader.vcf", skip = 1)

names(in_vcf) <- c("CHROM", "POS", "N_Alleles", "N_Chr", "REF", "ALT")
in_vcf$POS <- as.character(in_vcf$POS)
in_vcf$CHROM <- str_pad(in_vcf$CHROM, 2, pad = "0")
in_vcf$CHROM_POS <- apply(in_vcf[,c("CHROM", "POS")], 1, paste, collapse = "_")
in_vcf <- in_vcf[c("CHROM", "POS", "CHROM_POS", "N_Alleles", "N_Chr", "REF", "ALT")]
in_vcf$POS <- as.integer(in_vcf$POS)

REF <- in_frqs[in_frqs$REF == 1,]
REF <- REF[,c(1,2,3,seq(7,23,2))]

ALT <- in_frqs[in_frqs$ALT == 1,]
ALT <- ALT[,c(1,2,3,seq(6,22,2))]

names(REF) <- c("CHROM", "POS", "CHROM_POS", "Outgroup",
                "AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "D5", "A1")
names(ALT) <- c("CHROM", "POS", "CHROM_POS", "Outgroup",
                "AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "D5", "A1")

ancestral_vcf <- rbind(REF, ALT)

GERP_vcf <- ancestral_vcf[FALSE,]


for(i in c(01:26)){
  j <- str_pad(i, 2, pad = "0")
  print(j)
  GERP_name <- paste("../../../GERP_SCORES/roast", j, "msa.in.new.raw.rates", sep = ".")
  GERP <- read.table(GERP_name, header=FALSE)
  names(GERP) <- c('Neutral Rate','RS Score')
  temp_frqs <- ancestral_vcf[ancestral_vcf$CHROM == j,]
  temp_frqs$GERP <- GERP[temp_frqs$POS,2]
  GERP_vcf <- dplyr::bind_rows(temp_frqs, GERP_vcf)
}

print("GERP added")
BAD_results <- read.table("../all.GT.merged_results.txt", header = TRUE, sep = "\t")
BAD_results <- as.data.frame(BAD_results[,c(21, 15, 14, 3)])


print("BAD_results loaded")
GERP_vcf <- as.data.frame(GERP_vcf)
head(GERP_vcf)
print("")
head(BAD_results)
almost_vcf <- merge(x=GERP_vcf, y=BAD_results, by.x = "CHROM_POS", by.y = "SNP.Name", all.x = TRUE)

print("BAD_results and GERP merged")
synonymous <- read.table("../all.GT.se_variant.txt", header = FALSE)
names(synonymous) <- c("CHROM_POS", "Synonymous", "Syn_Gene")
final_vcf <- merge(x=almost_vcf, y=synonymous, by.x = "CHROM_POS", by.y = "CHROM_POS", all.x = TRUE)
print("merged BAD_GERP with Synonymous file")
write.table(final_vcf, file = "all.frqs.BAD.allGERP.pSONIC.noGC.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

GERP06 <- subset(final_vcf, (GERP > -2 & GERP != 0) | !is.na(GeneID) | !is.na(Synonymous))
write.table(GERP06, file = "GERP06_BAD.frqs.pSONIC.noGC.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

GERP_melt <- melt(GERP06, measure.vars = c("A1", "AD1", "AD3", "AD4", "AD5", "AD6",
                "AD7", "D5"), value.name = c("DAF"), variable.name = c("Population"))

#GERP_melt$SNP_ID_Pop <- apply(GERP_melt[,c("CHROM_POS", "Population")], 1, paste, collapse = "_")

write.table(GERP_melt, file="GERP06.melted.frqs.pSONIC.noGC.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

GERP_melt <- GERP_melt[c(GERP_melt$DAF > 0 & !is.na(GERP_melt$DAF)),]

write.table(GERP_melt, file="GERP06.melted.frqs.seg.pSONIC.noGC.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
