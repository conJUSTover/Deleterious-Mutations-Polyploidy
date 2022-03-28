library(stringr)
library(dplyr)
library(reshape2)

setwd("/work/LAS/jfw-lab/del_mut_gvcf/redo_f2/gVCFs/Australian/frqs/models/")

infile <- read.table("pSONIC.GT.NO_GC.GERP06_BAD.vcf", header = TRUE)

count_rows_subset <- function(col_1){
  return(2 * sum(c(col_1 == 2 & !is.na(col_1))))
}

Afile <- infile[c(infile$CHROM < 14 & !is.na(infile$CHROM)),]
Dfile <- infile[c(infile$CHROM > 13 & !is.na(infile$CHROM)),]

Arecessive <- data.frame(Population = names(Afile[,4:46]),
                         Mildly = apply(Afile[c(Afile$GERP <= 2 & Afile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(Afile[c(Afile$GERP <= 4 & Afile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(Afile[c(Afile$GERP <= 6 & Afile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(Afile[c(Afile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Diploid"),
                         Subgenome = rep("A"), NAs = rep("yes"))

Aadditive <- data.frame(Population = names(Afile[,4:46]),
                        Mildly = colSums(Afile[c(Afile$GERP <= 2 & Afile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(Afile[c(Afile$GERP <= 4 & Afile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(Afile[c(Afile$GERP <= 6 & Afile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(Afile[c(Afile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Diploid"),
                        Subgenome = rep("A"), NAs = rep("yes"))

Drecessive <- data.frame(Population = names(Dfile[,4:46]),
                         Mildly = apply(Dfile[c(Dfile$GERP <= 2 & Dfile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(Dfile[c(Dfile$GERP <= 4 & Dfile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(Dfile[c(Dfile$GERP <= 6 & Dfile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(Dfile[c(Dfile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Diploid"),
                         Subgenome = rep("D"), NAs = rep("yes"))

Dadditive <- data.frame(Population = names(Dfile[,4:46]),
                        Mildly = colSums(Dfile[c(Dfile$GERP <= 2 & Dfile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(Dfile[c(Dfile$GERP <= 4 & Dfile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(Dfile[c(Dfile$GERP <= 6 & Dfile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(Dfile[c(Dfile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Diploid"),
                        Subgenome = rep("D"), NAs = rep("yes"))

Dfile$averageDt <- rowMeans(Dfile[,c(4:41)], na.rm = TRUE)
DTfile <- Dfile[c(Dfile$averageDt < 1 & Dfile$averageDt > 0 & !is.na(Dfile$averageDt)),]

Afile$averageAt <- rowMeans(Afile[,c(4:37,42:46)], na.rm = TRUE)
ATfile <- Afile[c(Afile$averageAt < 1 & Afile$averageAt > 0 & !is.na(Afile$averageAt)),]

ATrecessive <- data.frame(Population = names(ATfile[,4:46]),
                         Mildly = apply(ATfile[c(ATfile$GERP <= 2 & ATfile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(ATfile[c(ATfile$GERP <= 4 & ATfile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(ATfile[c(ATfile$GERP <= 6 & ATfile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(ATfile[c(ATfile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Tetraploid"),
                         Subgenome = rep("A"), NAs = rep("yes"))

ATadditive <- data.frame(Population = names(ATfile[,4:46]),
                        Mildly = colSums(ATfile[c(ATfile$GERP <= 2 & ATfile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(ATfile[c(ATfile$GERP <= 4 & ATfile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(ATfile[c(ATfile$GERP <= 6 & ATfile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(ATfile[c(ATfile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Tetraploid"),
                        Subgenome = rep("A"), NAs = rep("yes"))

DTrecessive <- data.frame(Population = names(DTfile[,4:46]),
                         Mildly = apply(DTfile[c(DTfile$GERP <= 2 & DTfile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(DTfile[c(DTfile$GERP <= 4 & DTfile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(DTfile[c(DTfile$GERP <= 6 & DTfile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(DTfile[c(DTfile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Tetraploid"),
                         Subgenome = rep("D"), NAs = rep("yes"))

DTadditive <- data.frame(Population = names(DTfile[,4:46]),
                        Mildly = colSums(DTfile[c(DTfile$GERP <= 2 & DTfile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(DTfile[c(DTfile$GERP <= 4 & DTfile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(DTfile[c(DTfile$GERP <= 6 & DTfile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(DTfile[c(DTfile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Tetraploid"),
                        Subgenome = rep("D"), NAs = rep("yes"))

Afile <- Afile[complete.cases(Afile[c(4:37,42:46),]),]
Dfile <- Dfile[complete.cases(Dfile[c(4:41),]),]
ATfile <- ATfile[complete.cases(ATfile[c(4:37,42:46),]),]
DTfile <- DTfile[complete.cases(DTfile[c(4:41),]),]

ArecessiveNA <- data.frame(Population = names(Afile[,4:46]),
                         Mildly = apply(Afile[c(Afile$GERP <= 2 & Afile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(Afile[c(Afile$GERP <= 4 & Afile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(Afile[c(Afile$GERP <= 6 & Afile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(Afile[c(Afile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Diploid"),
                         Subgenome = rep("A"), NAs = rep("no"))

AadditiveNA <- data.frame(Population = names(Afile[,4:46]),
                        Mildly = colSums(Afile[c(Afile$GERP <= 2 & Afile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(Afile[c(Afile$GERP <= 4 & Afile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(Afile[c(Afile$GERP <= 6 & Afile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(Afile[c(Afile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Diploid"),
                        Subgenome = rep("A"), NAs = rep("no"))
 
DrecessiveNA <- data.frame(Population = names(Dfile[,4:46]),
                         Mildly = apply(Dfile[c(Dfile$GERP <= 2 & Dfile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(Dfile[c(Dfile$GERP <= 4 & Dfile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(Dfile[c(Dfile$GERP <= 6 & Dfile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(Dfile[c(Dfile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Diploid"),
                         Subgenome = rep("D"), NAs = rep("no"))

DadditiveNA <- data.frame(Population = names(Dfile[,4:46]),
                        Mildly = colSums(Dfile[c(Dfile$GERP <= 2 & Dfile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(Dfile[c(Dfile$GERP <= 4 & Dfile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(Dfile[c(Dfile$GERP <= 6 & Dfile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(Dfile[c(Dfile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Diploid"),
                        Subgenome = rep("D"), NAs = rep("no"))



ATrecessiveNA <- data.frame(Population = names(ATfile[,4:46]),
                         Mildly = apply(ATfile[c(ATfile$GERP <= 2 & ATfile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(ATfile[c(ATfile$GERP <= 4 & ATfile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(ATfile[c(ATfile$GERP <= 6 & ATfile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(ATfile[c(ATfile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Tetraploid"),
                         Subgenome = rep("A"), NAs = rep("no"))

ATadditiveNA <- data.frame(Population = names(ATfile[,4:46]),
                        Mildly = colSums(ATfile[c(ATfile$GERP <= 2 & ATfile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(ATfile[c(ATfile$GERP <= 4 & ATfile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(ATfile[c(ATfile$GERP <= 6 & ATfile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(ATfile[c(ATfile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Tetraploid"),
                        Subgenome = rep("A"), NAs = rep("no"))

DTrecessiveNA <- data.frame(Population = names(DTfile[,4:46]),
                         Mildly = apply(DTfile[c(DTfile$GERP <= 2 & DTfile$GERP > 0),4:46], 2, count_rows_subset),
                         Moderately = apply(DTfile[c(DTfile$GERP <= 4 & DTfile$GERP > 2),4:46], 2, count_rows_subset),
                         Highly = apply(DTfile[c(DTfile$GERP <= 6 & DTfile$GERP > 4),4:46], 2, count_rows_subset),
                         Synonymous = apply(DTfile[c(DTfile$Synonymous == "synonymous_variant"),4:46], 2, count_rows_subset),
                         Model = rep("Recessive"),
                         Filter = rep("Tetraploid"),
                         Subgenome = rep("D"), NAs = rep("no"))

DTadditiveNA <- data.frame(Population = names(DTfile[,4:46]),
                        Mildly = colSums(DTfile[c(DTfile$GERP <= 2 & DTfile$GERP > 0),4:46], na.rm = TRUE),
                        Moderately = colSums(DTfile[c(DTfile$GERP <= 4 & DTfile$GERP > 2),4:46], na.rm = TRUE),
                        Highly = colSums(DTfile[c(DTfile$GERP <= 6 & DTfile$GERP > 4),4:46], na.rm = TRUE),
                        Synonymous = colSums(DTfile[c(DTfile$Synonymous == "synonymous_variant"),4:46], na.rm = TRUE),
                        Model = rep("Additive"),
                        Filter = rep("Tetraploid"),
                        Subgenome = rep("D"), NAs = rep("no"))



models <- rbind(Arecessive, Aadditive, Drecessive, Dadditive, ATrecessive, ATadditive, DTrecessive, DTadditive, ArecessiveNA, AadditiveNA, DrecessiveNA, DadditiveNA, ATrecessiveNA, ATadditiveNA, DTrecessiveNA, DTadditiveNA)
write.table(models, file = "Add_Recessive_Models.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



