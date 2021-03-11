library(ggplot2)
library(reshape2)

setwd("~/GitHub/Deleterious-Mutations-Polyploidy/figures/sfs/")

Species <- c("A1", "AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "D5")


load_sfs <- function(pop, prefix){
  missense_AD1 <- t(read.table(paste(prefix,"missense.recode.", pop, ".sfs.txt", sep = ''), skip = 1))
  syn_AD1 <- t(read.table(paste(prefix,"synonymous.recode.", pop, ".sfs.txt", sep = ''), skip = 1))
  GERP_high_AD1 <- t(read.table(paste(prefix,"GERP_high_del.recode.", pop, ".sfs.txt", sep = ''), skip = 1))
  GERP_mod_AD1 <- t(read.table(paste(prefix,"GERP_mod_del.recode.", pop, ".sfs.txt", sep = ''), skip = 1))
  GERP_slight_AD1 <- t(read.table(paste(prefix,"GERP_slight_del.recode.", pop, ".sfs.txt", sep = ''), skip = 1))
  GERP_nn_AD1 <- t(read.table(paste(prefix,"GERP_nearly_neutral.recode.", pop, ".sfs.txt", sep = ''), skip = 1))
  BAD_AD1 <- t(read.table(paste(prefix,"BAD_deleterious.recode.", pop, ".sfs.txt", sep = ''), skip = 1))

  AD1 <- data.frame(Missense <- missense_AD1[,1],
                  Synonymous <- syn_AD1[,1],
                  BAD <- BAD_AD1[,1],
                  GERP_High <- GERP_high_AD1[,1],
                  GERP_mod <- GERP_mod_AD1[,1],
                  GERP_slight <- GERP_slight_AD1[,1],
                  GERP_nn <- GERP_nn_AD1[,1]) 
  names(AD1) <- c("Missense", "Synonymous", "BAD", "GERP_High", "GERP_Moderate", "GERP_Slightly", "GERP_Nearly_Neutral")

  AD1$DAF <- c(seq(0,1,by = 1/(nrow(AD1) -1)))
  AD1 <- AD1[c(AD1$DAF > 0 & AD1$DAF < 1),]

  AD1[,1:7] <- sweep(AD1[,1:7], 2, colSums(AD1[,1:7]), `/`)
  AD1_melt <- melt(AD1, id.vars = "DAF", variable.name = 'Type')
  AD1_melt$Plot <- c(rep("BAD", nrow(AD1)* 3), rep("GERP", nrow(AD1) * 4))
  return(AD1_melt)
}

make_plots <- function(species){
  trial <- load_sfs(species, "")
  trial$Genome <- rep("All")
  trial_At <- load_sfs(species, "At.")
  trial_At$Genome <- rep("At")
  trial_Dt <- load_sfs(species, "Dt.")
  trial_Dt$Genome <- rep("Dt")

  all <- rbind(trial, trial_At, trial_Dt)
  p <- ggplot(all, aes(fill = Type)) + 
    geom_col(aes(x = DAF, y = value), stat = "identity", position = "dodge") + 
    facet_grid(Plot~Genome) + theme_bw() + 
    theme(legend.position = "bottom", axis.title.y = element_blank(), 
          legend.title = element_blank(), legend.direction = "horizontal")
  return(p)
}


p <- make_plots("AD1")
png(paste("Figure_SFS_AD1.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

p <- make_plots("AD3")
png(paste("Figure_SFS_AD3.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

p <- make_plots("AD4")
png(paste("Figure_SFS_AD4.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

p <- make_plots("AD5")
png(paste("Figure_SFS_AD5.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

p <- make_plots("AD6")
png(paste("Figure_SFS_AD6.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

p <- make_plots("AD7")
png(paste("Figure_SFS_AD7.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

p <- make_plots("A1")
png(paste("Figure_SFS_A1.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

p <- make_plots("D5")
png(paste("Figure_SFS_D5.png", sep = ''), 6000, 2500, pointsize=12, res=600)
p
dev.off()

