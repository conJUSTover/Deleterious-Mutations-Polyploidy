library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
#library(patchwork)
setwd("~/Documents/GERP/good_filtering_results/")

nfile <- read.table("GERP06_BAD.frqs.pSONIC.NO_GC.txt", header = TRUE)
nfile <- nfile[complete.cases(nfile[,5:12]),]

homoeologs <- read.table("pSONIC.homoeologs.375_625.txt", header = TRUE)
homoeologs$At = substr(homoeologs$At,1,nchar(homoeologs$At)-2)
homoeologs$Dt = substr(homoeologs$Dt,1,nchar(homoeologs$Dt)-2)
homoeologs <- as.data.frame(sapply(homoeologs, function(v) {gsub("\\.","", v)}))

GERP_count_bygene <- function(df, sp){
  return(sum(df[,sp], na.rm = TRUE))
}

GERP_count_bygene_all <- function(gene, Species, low, high_inc){
  temp <- infile[c(infile$GERP > low & infile$GERP <= high_inc & !is.na(infile$GERP) & !is.na(infile$Syn_Gene) & infile$Syn_Gene == gene),]
  output <- sapply(Species, GERP_count_bygene, df=temp)
  return(output)
}

#At <- nfile[nfile$CHROM < 14,]
#Dt <- nfile[nfile$CHROM > 13,]
#At <- At[complete.cases(At[,c(5:10,12)]),]
#Dt <- Dt[complete.cases(Dt[,c(5:11)]),]
#At$average <- rowMeans(At[,c(5:10,12)], na.rm = TRUE)
#Dt$average <- rowMeans(Dt[,c(5:11)], na.rm = TRUE)
#At <- At[c(At$average > 0 & At$average < 1 & !is.na(At$average)),]
#Dt <- Dt[c(Dt$average > 0 & Dt$average < 1 & !is.na(Dt$average)),]
#infile <- rbind(At, Dt)
infile <- nfile
rm(nfile)
infile <- infile[complete.cases(infile[,c(5:10)]),]
infile$average <- rowMeans(infile[,c(5:10)], na.rm = TRUE)
infile <- infile[c(infile$average > 0 & infile$average < 1 & !is.na(infile$average)),]

rm(At)
rm(Dt)

print("Diploid Sites converted")

homoeologs[, c("At_AD1_Gerp02", "At_AD3_Gerp02", "At_AD4_Gerp02", "At_AD5_Gerp02", "At_AD6_Gerp02", "At_AD7_Gerp02", "At_A1_Gerp02", "At_D5_Gerp02")] <- t(apply(homoeologs, 1, function(x) GERP_count_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"), low=0, high_inc=2)))
homoeologs[, c("Dt_AD1_Gerp02", "Dt_AD3_Gerp02", "Dt_AD4_Gerp02", "Dt_AD5_Gerp02", "Dt_AD6_Gerp02", "Dt_AD7_Gerp02", "Dt_A1_Gerp02", "Dt_D5_Gerp02")] <- t(apply(homoeologs, 1, function(x) GERP_count_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"), low=0, high_inc=2)))
print("02 Sites completed")

homoeologs[, c("At_AD1_Gerp24", "At_AD3_Gerp24", "At_AD4_Gerp24", "At_AD5_Gerp24", "At_AD6_Gerp24", "At_AD7_Gerp24", "At_A1_Gerp24", "At_D5_Gerp24")] <- t(apply(homoeologs, 1, function(x) GERP_count_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"), low=2, high_inc=4)))
homoeologs[, c("Dt_AD1_Gerp24", "Dt_AD3_Gerp24", "Dt_AD4_Gerp24", "Dt_AD5_Gerp24", "Dt_AD6_Gerp24", "Dt_AD7_Gerp24", "Dt_A1_Gerp24", "Dt_D5_Gerp24")] <- t(apply(homoeologs, 1, function(x) GERP_count_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"), low=2, high_inc=4)))
print("24 Sites completed")

homoeologs[, c("At_AD1_Gerp46", "At_AD3_Gerp46", "At_AD4_Gerp46", "At_AD5_Gerp46", "At_AD6_Gerp46", "At_AD7_Gerp46", "At_A1_Gerp46", "At_D5_Gerp46")] <- t(apply(homoeologs, 1, function(x) GERP_count_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"), low=4, high_inc=6)))
homoeologs[, c("Dt_AD1_Gerp46", "Dt_AD3_Gerp46", "Dt_AD4_Gerp46", "Dt_AD5_Gerp46", "Dt_AD6_Gerp46", "Dt_AD7_Gerp46", "Dt_A1_Gerp46", "Dt_D5_Gerp46")] <- t(apply(homoeologs, 1, function(x) GERP_count_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"), low=4, high_inc=6)))

write.table(homoeologs, "Homoeologs_GERPCount_per_Pop.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
homoeologs <- read.table("Homoeologs_GERPCount_per_Pop.txt", header = TRUE)

new_sums <- data.frame(Sample = names(homoeologs[,c(3:50)]), Load = colSums(homoeologs[,c(3:50)]))
new_sums$Species  = rep(c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))
new_sums$Subgenome = rep(c("At", "Dt"), each = 8)
new_sums$Deleteriousness <- factor(rep(c("Mildly","Moderately","Strongly"), each = 16), levels = c("Mildly", "Moderately", "Strongly"))


new_sums <- new_sums[!c(new_sums$Species == "D5"),]
new_sums <- new_sums[!c(new_sums$Species == "A1"),]
#new_sums$Species  = rep(c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "Diploid"))

At <- new_sums[new_sums$Subgenome == "At",]
At$Species <- factor(At$Species, levels = c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7"))
Dt <- new_sums[new_sums$Subgenome == "Dt",]
Dt$Species <- factor(Dt$Species, levels = c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7"))
Dtold <- Dt[Dt$Species != "AD4",]
Dtnew <- Dt[Dt$Species == "AD4",]
Dt <- rbind(Dtold, Dtnew)
Atold <- At[At$Species != "AD4",]
Atnew <- At[At$Species == "AD4",]
At <- rbind(Atold, Atnew)

At2 <- At %>% 
  group_by(Deleteriousness) %>%
  mutate(New = Load/last(Load))

Dt2 <- Dt %>% 
  group_by(Deleteriousness) %>%
  mutate(New = Load/last(Load))

new_sums <- rbind(At2, Dt2)


ggplot(new_sums, aes(x = Deleteriousness, y = New, fill = Deleteriousness)) + 
  geom_bar(stat = "identity") + theme_bw() + scale_fill_grey(start=0.85, end=0) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "bottom") +
  facet_grid(Subgenome~Species, scales = "free_y") + ylab("Proportion of Diploid SNPs")

ggsave("Figure_GERP_Proportions_poly.pdf", width = 11, height = 6, units = "cm")
