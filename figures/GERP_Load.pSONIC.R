library(reshape2)
#library(ggplot2)
#library(patchwork)
setwd("~/Documents/GERP/good_filtering_results/")

infile <- read.table("GERP06_BAD.frqs.pSONIC.noGC.txt", header = TRUE)

homoeologs <- read.table("pSONIC.homoeologs.375_625.txt", header = TRUE)
homoeologs$At = substr(homoeologs$At,1,nchar(homoeologs$At)-2)
homoeologs$Dt = substr(homoeologs$Dt,1,nchar(homoeologs$Dt)-2)
homoeologs <- as.data.frame(sapply(homoeologs, function(v) {gsub("\\.","", v)}))

bonferonni <- 0.05 / 334097 #Number from other GERP_Load R script 

#count_BAD_gene <- function(gene){
#  return(nrow(infile[c(infile$Syn_Gene == gene & infile$MaskedP.value < bonferonni& !is.na(infile$Syn_Gene) & !is.na(infile$MaskedP.value)),]))
#}

#count_GERP <- function(gene){
#  return(nrow(infile[c(infile$Syn_Gene == gene & infile$GERP > 4 & !is.na(infile$Syn_Gene) & !is.na(infile$GERP)),]))
#}

GERP_load_At <- function(species){
  temp <- infile[c(infile$GERP > 0 & !is.na(infile$GERP) & !is.na(infile[,species]) & infile$CHROM < 14),]
  return(sum(temp[,species] * temp$GERP))
}

GERP_load_Dt <- function(species){
  temp <- infile[c(infile$GERP > 0 & !is.na(infile$GERP) & !is.na(infile[,species]) & infile$CHROM > 13),]
  return(sum(temp[,species] * temp$GERP))
}

GERP_load_genespace_At <- function(species){
  temp <- infile[c(infile$GERP > 0 & !is.na(infile$GERP) & !is.na(infile[,species]) & infile$CHROM < 14 & !is.na(infile$Syn_Gene)),]
  return(sum(temp[,species] * temp$GERP))
}

GERP_load_genespace_Dt <- function(species){
  temp <- infile[c(infile$GERP > 0 & !is.na(infile$GERP) & !is.na(infile[,species]) & infile$CHROM > 13 & !is.na(infile$Syn_Gene)),]
  return(sum(temp[,species] * temp$GERP))
}

BAD_load_genespace_At <- function(species){
  temp <- infile[c(infile$MaskedP.value < bonferonni & !is.na(infile$MaskedP.value) & !is.na(infile[,species]) & infile$CHROM < 14 & !is.na(infile$Syn_Gene)),]
  return(sum(temp[,species]))
}

BAD_load_genespace_Dt <- function(species){
  temp <- infile[c(infile$MaskedP.value < bonferonni & !is.na(infile$MaskedP.value) & !is.na(infile[,species]) & infile$CHROM > 13 & !is.na(infile$Syn_Gene)),]
  return(sum(temp[,species]))
}

GERP_bygene <- function(df, sp){
  return(sum(df[,sp] * df$GERP, na.rm = TRUE))
}

GERP_bygene_all <- function(gene, Species){
  temp <- infile[c(infile$GERP > 0 & !is.na(infile$GERP) & !is.na(infile$Syn_Gene) & infile$Syn_Gene == gene),]
  output <- sapply(Species, GERP_bygene, df=temp)
  return(output)
}

BAD_bygene <- function(df, sp){
  return(sum(df[,sp], na.rm = TRUE))
}

BAD_bygene_all <- function(gene, Species){
  temp <- infile[c(infile$MaskedP.value < bonferonni & !is.na(infile$MaskedP.value) & !is.na(infile$Syn_Gene) & infile$Syn_Gene == gene),]
  output <- sapply(Species, BAD_bygene, df=temp)
  return(output)
}


#GERP_loads <- data.frame(Species = c("AD1", "AD3", "AD4", 'AD5', 'AD6', 'AD7', "A1", 'D5'))

#GERP_loads$At_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
#GERP_loads$Dt_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
#GERP_loads$At_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
#GERP_loads$Dt_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
#GERP_loads$At_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
#GERP_loads$Dt_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

print("GERP Load done")


#homoeologs <- homoeologs[1:10,]


##GERP Load between homoeologs 
homoeologs[, c("At_AD1_Gerp_load", "At_AD3_Gerp_load", "At_AD4_Gerp_load", "At_AD5_Gerp_load", "At_AD6_Gerp_load", "At_AD7_Gerp_load", "At_A1_Gerp_load", "At_D5_Gerp_load")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_Gerp_load", "Dt_AD3_Gerp_load", "Dt_AD4_Gerp_load", "Dt_AD5_Gerp_load", "Dt_AD6_Gerp_load", "Dt_AD7_Gerp_load", "Dt_A1_Gerp_load", "Dt_D5_Gerp_load")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs[, c("At_AD1_BAD_load", "At_AD3_BAD_load", "At_AD4_BAD_load", "At_AD5_BAD_load", "At_AD6_BAD_load", "At_AD7_BAD_load", "At_A1_BAD_load", "At_D5_BAD_load")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_BAD_load", "Dt_AD3_BAD_load", "Dt_AD4_BAD_load", "Dt_AD5_BAD_load", "Dt_AD6_BAD_load", "Dt_AD7_BAD_load", "Dt_A1_BAD_load", "Dt_D5_BAD_load")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs BAD done")

####################################
###Convert to Diploid Seg Sites###
####################################

At <- infile[infile$CHROM < 14,]
Dt <- infile[infile$CHROM > 13,]
At$average <- rowMeans(At[,c(5:10,12)], na.rm = TRUE)
Dt$average <- rowMeans(Dt[,c(5:11)], na.rm = TRUE)
At <- At[c(At$average > 0 & At$average < 1 & !is.na(At$average)),]
Dt <- Dt[c(Dt$average > 0 & Dt$average < 1 & !is.na(Dt$average)),]
infile <- rbind(At, Dt)

rm(At)
rm(Dt)

print("Diploid Sites converted")

#infile$average <- rowMeans(infile[,6:11], na.rm = TRUE)
#infile <- infile[c(infile$average > 0 & infile$average < 1 & !is.na(infile$average)),]

#GERP_loads$At_diploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
#GERP_loads$Dt_diploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
#GERP_loads$At_diploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
#GERP_loads$Dt_diploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
#GERP_loads$At_diploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
#GERP_loads$Dt_diploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

#print("GERP Loads done")

##GERP Load between homoeologs 
homoeologs[, c("At_AD1_Gerp_load_dip", "At_AD3_Gerp_load_dip", "At_AD4_Gerp_load_dip", "At_AD5_Gerp_load_dip", "At_AD6_Gerp_load_dip", "At_AD7_Gerp_load_dip", "At_A1_Gerp_load_dip", "At_D5_Gerp_load_dip")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_Gerp_load_dip", "Dt_AD3_Gerp_load_dip", "Dt_AD4_Gerp_load_dip", "Dt_AD5_Gerp_load_dip", "Dt_AD6_Gerp_load_dip", "Dt_AD7_Gerp_load_dip", "Dt_A1_Gerp_load_dip", "Dt_D5_Gerp_load_dip")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs[, c("At_AD1_BAD_load_dip", "At_AD3_BAD_load_dip", "At_AD4_BAD_load_dip", "At_AD5_BAD_load_dip", "At_AD6_BAD_load_dip", "At_AD7_BAD_load_dip", "At_A1_BAD_load_dip", "At_D5_BAD_load_dip")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_BAD_load_dip", "Dt_AD3_BAD_load_dip", "Dt_AD4_BAD_load_dip", "Dt_AD5_BAD_load_dip", "Dt_AD6_BAD_load_dip", "Dt_AD7_BAD_load_dip", "Dt_A1_BAD_load_dip", "Dt_D5_BAD_load_dip")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs BAD done")

####################################
###Convert to Polyploid Seg Sites###
####################################

infile$average <- rowMeans(infile[,5:10], na.rm = TRUE)
infile <- infile[c(infile$average > 0 & infile$average < 1 & !is.na(infile$average)),]

#GERP_loads$At_polyploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
#GERP_loads$Dt_polyploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
#GERP_loads$At_polyploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
#GERP_loads$Dt_polyploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
#GERP_loads$At_polyploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
#GERP_loads$Dt_polyploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

print("Polyploid Sites converted")

##GERP Load between homoeologs 
homoeologs[, c("At_AD1_Gerp_load_poly", "At_AD3_Gerp_load_poly", "At_AD4_Gerp_load_poly", "At_AD5_Gerp_load_poly", "At_AD6_Gerp_load_poly", "At_AD7_Gerp_load_poly", "At_A1_Gerp_load_poly", "At_D5_Gerp_load_poly")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_Gerp_load_poly", "Dt_AD3_Gerp_load_poly", "Dt_AD4_Gerp_load_poly", "Dt_AD5_Gerp_load_poly", "Dt_AD6_Gerp_load_poly", "Dt_AD7_Gerp_load_poly", "Dt_A1_Gerp_load_poly", "Dt_D5_Gerp_load_poly")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs[, c("At_AD1_BAD_load_poly", "At_AD3_BAD_load_poly", "At_AD4_BAD_load_poly", "At_AD5_BAD_load_poly", "At_AD6_BAD_load_poly", "At_AD7_BAD_load_poly", "At_A1_BAD_load_poly", "At_D5_BAD_load_poly")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_BAD_load_poly", "Dt_AD3_BAD_load_poly", "Dt_AD4_BAD_load_poly", "Dt_AD5_BAD_load_poly", "Dt_AD6_BAD_load_poly", "Dt_AD7_BAD_load_poly", "Dt_A1_BAD_load_poly", "Dt_D5_BAD_load_poly")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs BAD done")


#####################################
###Convert to Population Seg Sites###
#####################################

infile[,5:12][infile[,5:12] == 1] <- 0

#GERP_loads$At_seg_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
#GERP_loads$Dt_seg_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
#GERP_loads$At_seg_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
#GERP_loads$Dt_seg_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
#GERP_loads$At_seg_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
#GERP_loads$Dt_segd_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

print("Segregating Sites converted")

##GERP Load between homoeologs 
homoeologs[, c("At_AD1_Gerp_load_seg", "At_AD3_Gerp_load_seg", "At_AD4_Gerp_load_seg", "At_AD5_Gerp_load_seg", "At_AD6_Gerp_load_seg", "At_AD7_Gerp_load_seg", "At_A1_Gerp_load_seg", "At_D5_Gerp_load_seg")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_Gerp_load_seg", "Dt_AD3_Gerp_load_seg", "Dt_AD4_Gerp_load_seg", "Dt_AD5_Gerp_load_seg", "Dt_AD6_Gerp_load_seg", "Dt_AD7_Gerp_load_seg", "Dt_A1_Gerp_load_seg", "Dt_D5_Gerp_load_seg")] <- t(apply(homoeologs, 1, function(x) GERP_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs[, c("At_AD1_BAD_load_seg", "At_AD3_BAD_load_seg", "At_AD4_BAD_load_seg", "At_AD5_BAD_load_seg", "At_AD6_BAD_load_seg", "At_AD7_BAD_load_seg", "At_A1_BAD_load_seg", "At_D5_BAD_load_seg")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[1], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
homoeologs[, c("Dt_AD1_BAD_load_seg", "Dt_AD3_BAD_load_seg", "Dt_AD4_BAD_load_seg", "Dt_AD5_BAD_load_seg", "Dt_AD6_BAD_load_seg", "Dt_AD7_BAD_load_seg", "Dt_A1_BAD_load_seg", "Dt_D5_BAD_load_seg")] <- t(apply(homoeologs, 1, function(x) BAD_bygene_all(x[2], c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))

print("Homoeologs BAD done")


##################
### Print Shit ###
##################

#write.table(GERP_loads, "Loads_Total_GERP_Bad.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(homoeologs, "Homoeologs_Load_per_Pop.pSNOIC.noGC.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
