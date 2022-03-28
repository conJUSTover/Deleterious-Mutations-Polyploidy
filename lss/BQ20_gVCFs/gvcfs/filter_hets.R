library(stringr)
library(reshape2)

setwd("/work/LAS/jfw-lab/del_mut_gvcf/BQ20_gVCFs/gvcfs")

#anc_pos <- read.table("all.Ancestral.POS.lst", header = TRUE)
#anc_pos$ChrPos <- apply(anc_pos[,c("CHROM", "POS")], 1, paste, collapse = "_")
#anc_pos <- anc_pos[anc_pos$CHROM == 1,]

hets_counts <- read.table(gzfile("all.hets.gz"), skip = 1)
hets_counts$V1 <- str_pad(hets_counts$V1, 2, pad = "0")
#hets_counts$ChrPos <- apply(hets_counts[,c("V1", "V2")], 1, paste, collapse = "_")



hets_per <- data.frame(Chr = as.integer(hets_counts$V1), POS = hets_counts$V2, 
                       ChrPos = apply(hets_counts[,c("V1", "V2")], 1, paste, collapse = "_"))

hets_per$AD1W <- hets_counts$V3 / hets_counts$V4
hets_per$AD1L1 <- hets_counts$V5 / hets_counts$V6
hets_per$AD1L2 <- hets_counts$V7 / hets_counts$V8
hets_per$AD1D <- hets_counts$V9 / hets_counts$V10
hets_per$AD2L1 <- hets_counts$V11 / hets_counts$V12
hets_per$AD2L2 <- hets_counts$V13 / hets_counts$V14
hets_per$AD2D <- hets_counts$V15 / hets_counts$V16
hets_per$AD3 <- hets_counts$V17 / hets_counts$V18
hets_per$AD4 <- hets_counts$V19 / hets_counts$V20
hets_per$AD5 <- hets_counts$V21 / hets_counts$V22
hets_per$AD6 <- hets_counts$V23 / hets_counts$V24
hets_per$AD7 <- hets_counts$V25 / hets_counts$V26
hets_per$A1W <- hets_counts$V27 / hets_counts$V28
hets_per$A1D <- hets_counts$V29 / hets_counts$V30
hets_per$A2W <- hets_counts$V31 / hets_counts$V32
hets_per$D5 <- hets_counts$V33 / hets_counts$V34

rm(hets_counts)

hets_melt <- hets_per[hets_per$Chr == 1,]
hets_melt <- melt(hets_melt, measure.vars = c("A1D", "A1W", "A2W", "AD1D",
                "AD1L1", "AD1L2", "AD1W", "AD2D", "AD2L1", "AD2L2", "AD3",
                "AD4", "AD5", "AD6", "AD7", "D5"), value.name = c("Hets"),
                variable.name = c("Population"))
hets_melt <- hets_melt[c(hets_melt$Hets == 1),]


for(i in 2:26){
  chr_per <- hets_per[hets_per$Chr == i,]
  chr_per <- melt(chr_per, measure.vars = c("A1D", "A1W", "A2W", "AD1D", 
                "AD1L1", "AD1L2", "AD1W", "AD2D", "AD2L1", "AD2L2", "AD3", 
                "AD4", "AD5", "AD6", "AD7", "D5"), value.name = c("Hets"), 
                variable.name = c("Population"))
  chr_per <- chr_per[c(chr_per$Hets == 1),]
  hets_melt <- dplyr::bind_rows(chr_per, hets_melt)
  rm(chr_per)
}

hets_melt$CHR_POS_POS <- apply(hets_melt[,c("ChrPos", "Population")], 1, paste, collapse = "_")

write.table(hets_melt, file="All_Hets.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

