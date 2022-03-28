library(stringr)
library(dplyr)
library(reshape2)

setwd("/work/LAS/jfw-lab/del_mut_gvcf/BQ20_gVCFs/gvcfs/new_filtering/frqs")

infile <- read.table("all.GT.GERP06_BAD.vcf", header = TRUE)

count_rows_subset <- function(col_1){
  return(2 * sum(c(col_1 == 2 & !is.na(col_1))))
}

Afile <- infile[c(infile$CHROM < 14 & !is.na(infile$CHROM)),]
Dfile <- infile[c(infile$CHROM > 13 & !is.na(infile$CHROM)),]

Arecessive <- data.frame(Population = names(Afile[,4:128]),
                         Mildly = apply(Afile[c(Afile$GERP <= 2 & Afile$GERP > 0),4:128], 2, count_rows_subset),
                         Moderately = apply(Afile[c(Afile$GERP <= 4 & Afile$GERP > 2),4:128], 2, count_rows_subset),
                         Highly = apply(Afile[c(Afile$GERP <= 6 & Afile$GERP > 4),4:128], 2, count_rows_subset),
                         Synonymous = apply(Afile[c(Afile$Synonymous == "synonymous_variant"),4:128], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Diploid"),
                         Subgenome = rep("A"))

Aadditive <- data.frame(Population = names(Afile[,4:128]),
                        Mildly = colSums(Afile[c(Afile$GERP <= 2 & Afile$GERP > 0),4:128], na.rm = TRUE),
                        Moderately = colSums(Afile[c(Afile$GERP <= 4 & Afile$GERP > 2),4:128], na.rm = TRUE),
                        Highly = colSums(Afile[c(Afile$GERP <= 6 & Afile$GERP > 4),4:128], na.rm = TRUE),
                        Synonymous = colSums(Afile[c(Afile$Synonymous == "synonymous_variant"),4:128], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Diploid"),
                        Subgenome = rep("A"))

Drecessive <- data.frame(Population = names(Dfile[,4:128]),
                         Mildly = apply(Dfile[c(Dfile$GERP <= 2 & Dfile$GERP > 0),4:128], 2, count_rows_subset),
                         Moderately = apply(Dfile[c(Dfile$GERP <= 4 & Dfile$GERP > 2),4:128], 2, count_rows_subset),
                         Highly = apply(Dfile[c(Dfile$GERP <= 6 & Dfile$GERP > 4),4:128], 2, count_rows_subset),
                         Synonymous = apply(Dfile[c(Dfile$Synonymous == "synonymous_variant"),4:128], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Diploid"),
                         Subgenome = rep("D"))

Dadditive <- data.frame(Population = names(Dfile[,4:128]),
                        Mildly = colSums(Dfile[c(Dfile$GERP <= 2 & Dfile$GERP > 0),4:128], na.rm = TRUE),
                        Moderately = colSums(Dfile[c(Dfile$GERP <= 4 & Dfile$GERP > 2),4:128], na.rm = TRUE),
                        Highly = colSums(Dfile[c(Dfile$GERP <= 6 & Dfile$GERP > 4),4:128], na.rm = TRUE),
                        Synonymous = colSums(Dfile[c(Dfile$Synonymous == "synonymous_variant"),4:128], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Diploid"),
                        Subgenome = rep("D"))

Dfile$averageDt <- rowMeans(Dfile[,c(4:109,114:128)], na.rm = TRUE)
DTfile <- Dfile[c(Dfile$averageDt < 1 & Dfile$averageDt > 0 & !is.na(Dfile$averageDt)),]

Afile$averageAt <- rowMeans(Afile[,c(4:5,16:113)], na.rm = TRUE)
ATfile <- Afile[c(Afile$averageAt < 1 & Afile$averageAt > 0 & !is.na(Afile$averageAt)),]

ATrecessive <- data.frame(Population = names(ATfile[,4:128]),
                         Mildly = apply(ATfile[c(ATfile$GERP <= 2 & ATfile$GERP > 0),4:128], 2, count_rows_subset),
                         Moderately = apply(ATfile[c(ATfile$GERP <= 4 & ATfile$GERP > 2),4:128], 2, count_rows_subset),
                         Highly = apply(ATfile[c(ATfile$GERP <= 6 & ATfile$GERP > 4),4:128], 2, count_rows_subset),
                         Synonymous = apply(ATfile[c(ATfile$Synonymous == "synonymous_variant"),4:128], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Tetraploid"),
                         Subgenome = rep("A"))

ATadditive <- data.frame(Population = names(ATfile[,4:128]),
                        Mildly = colSums(ATfile[c(ATfile$GERP <= 2 & ATfile$GERP > 0),4:128], na.rm = TRUE),
                        Moderately = colSums(ATfile[c(ATfile$GERP <= 4 & ATfile$GERP > 2),4:128], na.rm = TRUE),
                        Highly = colSums(ATfile[c(ATfile$GERP <= 6 & ATfile$GERP > 4),4:128], na.rm = TRUE),
                        Synonymous = colSums(ATfile[c(ATfile$Synonymous == "synonymous_variant"),4:128], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Tetraploid"),
                        Subgenome = rep("A"))

DTrecessive <- data.frame(Population = names(DTfile[,4:128]),
                         Mildly = apply(DTfile[c(DTfile$GERP <= 2 & DTfile$GERP > 0),4:128], 2, count_rows_subset),
                         Moderately = apply(DTfile[c(DTfile$GERP <= 4 & DTfile$GERP > 2),4:128], 2, count_rows_subset),
                         Highly = apply(DTfile[c(DTfile$GERP <= 6 & DTfile$GERP > 4),4:128], 2, count_rows_subset),
                         Synonymous = apply(DTfile[c(DTfile$Synonymous == "synonymous_variant"),4:128], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Tetraploid"),
                         Subgenome = rep("D"))

DTadditive <- data.frame(Population = names(DTfile[,4:128]),
                        Mildly = colSums(DTfile[c(DTfile$GERP <= 2 & DTfile$GERP > 0),4:128], na.rm = TRUE),
                        Moderately = colSums(DTfile[c(DTfile$GERP <= 4 & DTfile$GERP > 2),4:128], na.rm = TRUE),
                        Highly = colSums(DTfile[c(DTfile$GERP <= 6 & DTfile$GERP > 4),4:128], na.rm = TRUE),
                        Synonymous = colSums(DTfile[c(DTfile$Synonymous == "synonymous_variant"),4:128], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Tetraploid"),
                        Subgenome = rep("D"))



models <- rbind(Arecessive, Aadditive, Drecessive, Dadditive, ATrecessive, ATadditive, DTrecessive, DTadditive)
write.table(models, file = "Add_Recessive_Models.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



