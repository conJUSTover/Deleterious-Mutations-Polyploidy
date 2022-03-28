#library(reshape2)
#library(ggplot2)
#library(patchwork)
setwd("/work/LAS/jfw-lab/del_mut_gvcf/redo_f2/gVCFs/Australian/frqs/")

infile <- read.table("GERP06_BAD.frqs.txt", header = TRUE)

homoeologs <- read.table("pSONIC.homoeologs.375_625.txt", header = TRUE)
homoeologs$At = substr(homoeologs$At,1,nchar(homoeologs$At)-2)
homoeologs$Dt = substr(homoeologs$Dt,1,nchar(homoeologs$Dt)-2)
homoeologs <- as.data.frame(sapply(homoeologs, function(v) {gsub("\\.","", v)}))

bonferonni <- 0.05 / nrow(infile[!is.na(infile$MaskedConstraint),])

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

GERP_load_bygene <- function(gene, species){
  temp <- infile[c(infile$GERP > 0 & !is.na(infile$GERP) & !is.na(infile$Syn_Gene) & infile$Syn_Gene == gene),]
  return(sum(temp[,species] * temp$GERP))
}

BAD_load_bygene <- function(gene, species){
  temp <- infile[c(infile$MaskedP.value < bonferonni & !is.na(infile$MaskedP.value) & !is.na(infile$Syn_Gene) & infile$Syn_Gene == gene),]
  return(sum(temp[,species]))
}


GERP_loads <- data.frame(Species = c("AD1", "AD3", "AD4", 'AD5', 'AD6', 'AD7', "A1", 'D5'))

GERP_loads$At_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
GERP_loads$Dt_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
GERP_loads$At_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
GERP_loads$Dt_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
GERP_loads$At_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
GERP_loads$Dt_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

print("GERP Load done")
#homoeologs <- homoeologs[1:10,]

#homoeologs$At_BAD_muts <- unlist(lapply(homoeologs$At, count_BAD_gene))
#homoeologs$Dt_BAD_muts <- unlist(lapply(homoeologs$Dt, count_BAD_gene))
#homoeologs$Dt_At_diff <- homoeologs$Dt_BAD_muts - homoeologs$At_BAD_muts
#sum(homoeologs$At_BAD_muts)
#sum(homoeologs$Dt_BAD_muts)

#homoeologs$At_GERP_High <- unlist(lapply(homoeologs$At, count_GERP))
#homoeologs$Dt_GERP_High <- unlist(lapply(homoeologs$Dt, count_GERP))
#homoeologs$Dt_At_GERP_High_diff <- homoeologs$Dt_GERP_High - homoeologs$At_GERP_High


##GERP Load between homoeologs 
homoeologs$At_AD1_Gerp_load <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_Gerp_load <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD1"))
homoeologs$At_AD3_Gerp_load <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_Gerp_load <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD3"))
homoeologs$At_AD4_Gerp_load <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_Gerp_load <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD4"))
homoeologs$At_AD5_Gerp_load <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_Gerp_load <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD5"))
homoeologs$At_AD6_Gerp_load <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_Gerp_load <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD6"))
homoeologs$At_AD7_Gerp_load <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_Gerp_load <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD7"))
homoeologs$At_A1_Gerp_load <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "A1"))
homoeologs$Dt_D5_Gerp_load <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "D5"))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs$At_AD1_BAD_load <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_BAD_load <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD1"))
homoeologs$At_AD3_BAD_load <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_BAD_load <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD3"))
homoeologs$At_AD4_BAD_load <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_BAD_load <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD4"))
homoeologs$At_AD5_BAD_load <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_BAD_load <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD5"))
homoeologs$At_AD6_BAD_load <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_BAD_load <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD6"))
homoeologs$At_AD7_BAD_load <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_BAD_load <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD7"))
homoeologs$At_A1_BAD_load <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "A1"))
homoeologs$Dt_D5_BAD_load <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "D5"))

print("Homoeologs BAD done")

####################################
###Convert to Diploid Seg Sites###
####################################

At <- infile[infile$CHROM < 14,]
Dt <- infile[infile$CHROM > 13,]
At$average <- rowMeans(At[,c(5:10,12)], na.rm = TRUE)
Dt$average <- rowMeans(Dt[,c(5:11)], na.rm = TRUE)
At <- At[c(At$average > 0 & At$average < 1 & !is.na(At$average)),]
Dt <- At[c(Dt$average > 0 & Dt$average < 1 & !is.na(Dt$average)),]
infile <- rbind(At, Dt)

print("Diploid Sites converted")

#infile$average <- rowMeans(infile[,6:11], na.rm = TRUE)
#infile <- infile[c(infile$average > 0 & infile$average < 1 & !is.na(infile$average)),]

GERP_loads$At_diploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
GERP_loads$Dt_diploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
GERP_loads$At_diploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
GERP_loads$Dt_diploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
GERP_loads$At_diploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
GERP_loads$Dt_diploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

print("GERP Loads done")

##GERP Load between homoeologs 
homoeologs$At_AD1_Gerp_load_dip <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_Gerp_load_dip <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD1"))
homoeologs$At_AD3_Gerp_load_dip <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_Gerp_load_dip <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD3"))
homoeologs$At_AD4_Gerp_load_dip <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_Gerp_load_dip <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD4"))
homoeologs$At_AD5_Gerp_load_dip <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_Gerp_load_dip <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD5"))
homoeologs$At_AD6_Gerp_load_dip <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_Gerp_load_dip <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD6"))
homoeologs$At_AD7_Gerp_load_dip <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_Gerp_load_dip <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD7"))
homoeologs$At_A1_Gerp_load_dip <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "A1"))
homoeologs$Dt_D5_Gerp_load_dip <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "D5"))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs$At_AD1_BAD_load_dip <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_BAD_load_dip <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD1"))
homoeologs$At_AD3_BAD_load_dip <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_BAD_load_dip <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD3"))
homoeologs$At_AD4_BAD_load_dip <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_BAD_load_dip <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD4"))
homoeologs$At_AD5_BAD_load_dip <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_BAD_load_dip <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD5"))
homoeologs$At_AD6_BAD_load_dip <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_BAD_load_dip <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD6"))
homoeologs$At_AD7_BAD_load_dip <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_BAD_load_dip <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD7"))
homoeologs$At_A1_BAD_load_dip <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "A1"))
homoeologs$Dt_D5_BAD_load_dip <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "D5"))

print("Homoeologs BAD done")

####################################
###Convert to Polyploid Seg Sites###
####################################

infile$average <- rowMeans(infile[,5:10], na.rm = TRUE)
infile <- infile[c(infile$average > 0 & infile$average < 1 & !is.na(infile$average)),]

GERP_loads$At_polyploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
GERP_loads$Dt_polyploid_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
GERP_loads$At_polyploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
GERP_loads$Dt_polyploid_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
GERP_loads$At_polyploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
GERP_loads$Dt_polyploid_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

print("Polyploid Sites converted")

##GERP Load between homoeologs 
homoeologs$At_AD1_Gerp_load_poly <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_Gerp_load_poly <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD1"))
homoeologs$At_AD3_Gerp_load_poly <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_Gerp_load_poly <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD3"))
homoeologs$At_AD4_Gerp_load_poly <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_Gerp_load_poly <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD4"))
homoeologs$At_AD5_Gerp_load_poly <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_Gerp_load_poly <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD5"))
homoeologs$At_AD6_Gerp_load_poly <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_Gerp_load_poly <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD6"))
homoeologs$At_AD7_Gerp_load_poly <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_Gerp_load_poly <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD7"))
homoeologs$At_A1_Gerp_load_poly <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "A1"))
homoeologs$Dt_D5_Gerp_load_poly <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "D5"))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs$At_AD1_BAD_load_poly <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_BAD_load_poly <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD1"))
homoeologs$At_AD3_BAD_load_poly <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_BAD_load_poly <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD3"))
homoeologs$At_AD4_BAD_load_poly <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_BAD_load_poly <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD4"))
homoeologs$At_AD5_BAD_load_poly <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_BAD_load_poly <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD5"))
homoeologs$At_AD6_BAD_load_poly <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_BAD_load_poly <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD6"))
homoeologs$At_AD7_BAD_load_poly <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_BAD_load_poly <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD7"))
homoeologs$At_A1_BAD_load_poly <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "A1"))
homoeologs$Dt_D5_BAD_load_poly <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "D5"))

print("Homoeologs BAD done")

#####################################
###Convert to Population Seg Sites###
#####################################

infile[,5:12][infile[,5:12] == 1] <- 0

GERP_loads$At_seg_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_At))
GERP_loads$Dt_seg_model_genome <- unlist(lapply(GERP_loads$Species, GERP_load_Dt))
GERP_loads$At_seg_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_At))
GERP_loads$Dt_seg_model_genespace <- unlist(lapply(GERP_loads$Species, GERP_load_genespace_Dt))
GERP_loads$At_seg_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_At))
GERP_loads$Dt_segd_model_genespace_BAD <- unlist(lapply(GERP_loads$Species, BAD_load_genespace_Dt))

print("Segregating Sites converted")

##GERP Load between homoeologs 
homoeologs$At_AD1_Gerp_load_seg <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_Gerp_load_seg <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD1"))
homoeologs$At_AD3_Gerp_load_seg <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_Gerp_load_seg <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD3"))
homoeologs$At_AD4_Gerp_load_seg <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_Gerp_load_seg <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD4"))
homoeologs$At_AD5_Gerp_load_seg <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_Gerp_load_seg <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD5"))
homoeologs$At_AD6_Gerp_load_seg <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_Gerp_load_seg <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD6"))
homoeologs$At_AD7_Gerp_load_seg <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_Gerp_load_seg <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "AD7"))
homoeologs$At_A1_Gerp_load_seg <- unlist(lapply(homoeologs$At, GERP_load_bygene, species = "A1"))
homoeologs$Dt_D5_Gerp_load_seg <- unlist(lapply(homoeologs$Dt, GERP_load_bygene, species = "D5"))

print("Homoeologs GERP done")

###BAD Load between homoeologs 
homoeologs$At_AD1_BAD_load_seg <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD1"))
homoeologs$Dt_AD1_BAD_load_seg <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD1"))
homoeologs$At_AD3_BAD_load_seg <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD3"))
homoeologs$Dt_AD3_BAD_load_seg <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD3"))
homoeologs$At_AD4_BAD_load_seg <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD4"))
homoeologs$Dt_AD4_BAD_load_seg <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD4"))
homoeologs$At_AD5_BAD_load_seg <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD5"))
homoeologs$Dt_AD5_BAD_load_seg <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD5"))
homoeologs$At_AD6_BAD_load_seg <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD6"))
homoeologs$Dt_AD6_BAD_load_seg <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD6"))
homoeologs$At_AD7_BAD_load_seg <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "AD7"))
homoeologs$Dt_AD7_BAD_load_seg <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "AD7"))
homoeologs$At_A1_BAD_load_seg <- unlist(lapply(homoeologs$At, BAD_load_bygene, species = "A1"))
homoeologs$Dt_D5_BAD_load_seg <- unlist(lapply(homoeologs$Dt, BAD_load_bygene, species = "D5"))

print("Homoeologs BAD done")

##################
### Print Shit ###
##################

write.table(GERP_loads, "Loads_Total_GERP_Bad.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(homoeologs, "Homoeologs_Load_per_Pop.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
