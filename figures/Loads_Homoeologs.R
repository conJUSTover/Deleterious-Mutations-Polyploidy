library(ggplot2)
library(ggpubr)
library(reshape2)
library(patchwork)
library(RColorBrewer)

setwd("~/Documents/GERP/good_filtering_results/")

homoeolog_loads <- read.table("Homoeologs_Load_per_Pop.txt", header = TRUE)
load_sums <- data.frame(Sample = names(homoeolog_loads[,c(3:258)]), Load = colSums(homoeolog_loads[,c(3:258)]))
load_sums$Species  = rep(c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))
load_sums$Phylogeny = rep(c("Basal", "Diploid", "Polyploid", "Seg"), each = 64)
load_sums$Score = rep(c(rep("GERP", times = 16), rep("BAD", times = 16), rep("Synonymous", times=16), rep("Nonsynonymous", times=16)))
load_sums$Subgenome = rep(c("At", "Dt"), each = 8)
load_sums$Scale = rep("Homoeologs")

BAD <- load_sums[load_sums$Score == "BAD",]
Nonsynon <- load_sums[load_sums$Score == "Nonsynonymous",]
BAD$percent <- (BAD$Load / Nonsynon$Load) * 100
BAD$Load <- BAD$percent
BAD <- BAD[,-c(8)]
BAD$Score <- rep("Percent BAD /\nNonsynonymous")

homoeolog_loads <- read.table("Homoeologs_Load_per_Pop.pSNOIC.NO_GC.txt", header = TRUE)
new_sums <- data.frame(Sample = names(homoeolog_loads[,c(3:130)]), Load = colSums(homoeolog_loads[,c(3:130)]))
new_sums$Species  = rep(c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))
new_sums$Phylogeny = rep(c("Basal", "Diploid", "Polyploid", "Seg"), each = 32)
new_sums$Score = rep(c(rep("GERP", times = 16), rep("BAD", times = 16)))
new_sums$Subgenome = rep(c("At", "Dt"), each = 8)
new_sums$Scale = rep("Homoeologs\n(no GC)")

#output <- rbind(load_sums, new_sums)
output <- rbind(load_sums, BAD)
output <- output[!c(output$Species == "A1" & output$Subgenome == "Dt"),]
output <- output[!c(output$Species == "D5" & output$Subgenome == "At"),]
output <- output[!c(output$Species == "A1" & output$Phylogeny == "Polyploid"),]
output <- output[!c(output$Species == "D5" & output$Phylogeny == "Polyploid"),]



P <- output[c(output$Score == "GERP"),] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 2) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_brewer(palette = "Accent")+ ylim(0.0, NA)
P$layers <- c(geom_boxplot(outlier.shape = NA), P$layers)

P <- P + theme(legend.position = "none", axis.title.x = element_blank(), strip.text.y = element_blank()) #+ 
  ylab("GERP Load") 
#P

Q <- output[output$Score == "BAD",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 2) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_brewer(palette = "Accent") + ylim(0.0, NA) + 
  scale_y_continuous(position = "right") 
Q$layers <- c(geom_boxplot(outlier.shape = NA), Q$layers)

Q <- Q + theme(axis.title.x = element_blank(), strip.text.y = element_blank()) #+ 
  ylab("BAD_Mutations Load") 
#Q

R <- output[output$Score == "Synonymous",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 2) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_brewer(palette = "Accent")  + ylim(0.0, NA) 
R$layers <- c(geom_boxplot(outlier.shape = NA), R$layers)
R <- R + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank()) 

S <- output[output$Score == "Nonsynonymous",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 2) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_brewer(palette = "Accent") + ylim(0.0, NA) 
S$layers <- c(geom_boxplot(outlier.shape = NA), S$layers)
S <- S + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank())  

U <- output[output$Score == "Percent BAD /\nNonsynonymous",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 2) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_brewer(palette = "Accent") #+ ylim(0.0, NA) 
U$layers <- c(geom_boxplot(outlier.shape = NA), U$layers)
U <- U + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank())  

png("Figure_Loads_per_Phylogeny_percentages.Homoeologs.png", pointsize = 12, width = 5000, height = 5000, res = 600)
U + R + S + P + Q + plot_layout(widths = c(1,1,1,1,1))
dev.off()
