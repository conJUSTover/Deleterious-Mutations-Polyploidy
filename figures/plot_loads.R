library(ggplot2)
library(ggpubr)
library(reshape2)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(zoo)

setwd("~/Documents/GERP/good_filtering_results/")
infile <- read.table("Loads_Total_GERP_Bad.txt", header = TRUE)

my.cols <- brewer.pal(9, "Set1")
my.cols[9] <- "#000000"
my.cols <- my.cols[c(1,2,4:9)]
#homoeolog_loads <- read.table("Homoeologs_Load_per_Pop.txt", header = TRUE)

genome <- infile[,c(1:3)]
names(genome) <- c("Species", "At", "Dt")
output <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
output$Scale <- rep("Genome")
output$Phylogeny <- rep("Basal")
#output <- output[!c(output$Species == "A1" & output$Subgenome == "Dt"),]
#output <- output[!c(output$Species == "D5" & output$Subgenome == "At"),]
output$Score <- rep("GERP\nLoad")

genome <- infile[,c(1,4:5)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Basal")
temp$Score <- rep("GERP\nLoad")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)

genome <- infile[,c(1,6:7)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Basal")
temp$Score <- rep("BAD_Mutations\nLoad")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)

genome <- infile[,c(1,8:9)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Basal")
temp$Score <- rep("Synonymous\nMutations")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)


genome <- infile[,c(1,10:11)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Basal")
temp$Score <- rep("Nonsynonymous\nMutations")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)

##Diploids
genome <- infile[,c(1,12:13)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genome")
temp$Phylogeny <- rep("Diploid")
temp$Score <- rep("GERP\nLoad")
temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)

genome <- infile[,c(1,14:15)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Diploid")
temp$Score <- rep("GERP\nLoad")
temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)

genome <- infile[,c(1,16:17)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Diploid")
temp$Score <- rep("BAD_Mutations\nLoad")
temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)

genome <- infile[,c(1,18:19)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Diploid")
temp$Score <- rep("Synonymous\nMutations")
temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)

genome <- infile[,c(1,20:21)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Diploid")
temp$Score <- rep("Nonsynonymous\nMutations")
temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
output <- rbind(output, temp)


##Polyploid 
genome <- infile[,c(1,22:23)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genome")
temp$Phylogeny <- rep("Polyploid")
temp$Score <- rep("GERP\nLoad")
temp <- temp[!c(temp$Species == "D5"),]
temp <- temp[!c(temp$Species == "A1"),]
output <- rbind(output, temp)

genome <- infile[,c(1,24:25)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Polyploid")
temp$Score <- rep("GERP\nLoad")
temp <- temp[!c(temp$Species == "D5"),]
temp <- temp[!c(temp$Species == "A1"),]
output <- rbind(output, temp)

genome <- infile[,c(1,26:27)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Polyploid")
temp$Score <- rep("BAD_Mutations\nLoad")
temp <- temp[!c(temp$Species == "D5"),]
temp <- temp[!c(temp$Species == "A1"),]
output <- rbind(output, temp)

genome <- infile[,c(1,28:29)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Polyploid")
temp$Score <- rep("Synonymous\nMutations")
temp <- temp[!c(temp$Species == "A1"),]
temp <- temp[!c(temp$Species == "D5"),]
output <- rbind(output, temp)

genome <- infile[,c(1,30:31)]
names(genome) <- c("Species", "At", "Dt")
temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
temp$Scale <- rep("Genespace")
temp$Phylogeny <- rep("Polyploid")
temp$Score <- rep("Nonsynonymous\nMutations")
temp <- temp[!c(temp$Species == "A1"),]
temp <- temp[!c(temp$Species == "D5"),]
output <- rbind(output, temp)


#Segregating
#genome <- infile[,c(1,32:33)]
#names(genome) <- c("Species", "At", "Dt")
#temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
#temp$Scale <- rep("Genome")
#temp$Phylogeny <- rep("Segregating")
#temp$Score <- rep("GERP\nLoad")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
#output <- rbind(output, temp)

#genome <- infile[,c(1,34:35)]
#names(genome) <- c("Species", "At", "Dt")
#temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
#temp$Scale <- rep("Genespace")
#temp$Phylogeny <- rep("Segregating")
#temp$Score <- rep("GERP\nLoad")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
#output <- rbind(output, temp)

#genome <- infile[,c(1,36:37)]
#names(genome) <- c("Species", "At", "Dt")
#temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
#temp$Scale <- rep("Genespace")
#temp$Phylogeny <- rep("Segregating")
#temp$Score <- rep("BAD_Mutations\nLoad")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
#output <- rbind(output, temp)

#genome <- infile[,c(1,38:39)]
#names(genome) <- c("Species", "At", "Dt")
#temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
#temp$Scale <- rep("Genespace")
#temp$Phylogeny <- rep("Segregating")
#temp$Score <- rep("Synonymous\nMutations")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
#output <- rbind(output, temp)


#genome <- infile[,c(1,40:41)]
#names(genome) <- c("Species", "At", "Dt")
#temp <- melt(genome, measure.vars = c("At", "Dt"), variable.name = "Subgenome", value.name = "Load")
#temp$Scale <- rep("Genespace")
#temp$Phylogeny <- rep("Segregating")
#temp$Score <- rep("Nonsynonymous\nMutations")
#temp <- temp[!c(temp$Species == "A1" & temp$Subgenome == "Dt"),]
#temp <- temp[!c(temp$Species == "D5" & temp$Subgenome == "At"),]
#output <- rbind(output, temp)

output$Scale <- factor(output$Scale, levels = c("Genome", "Genespace"), labels = c("Genome", "Genespace"))
output <- output[output$Scale == "Genespace",]
output$Score <- factor(output$Score, levels = c("Synonymous\nMutations", "Nonsynonymous\nMutations", "GERP\nLoad", "BAD_Mutations\nLoad"))


BAD <- output[output$Score == "BAD_Mutations\nLoad",]
Nonsynon <- output[output$Score == "Nonsynonymous\nMutations",]
BAD$percent <- (BAD$Load / Nonsynon$Load) * 100
BAD$Load <- BAD$percent
BAD <- BAD[,-c(7)]
BAD$Score <- rep("All Genes")
Percent_Bad <- BAD
###Plot Everything 

P <- output[c(output$Score == "GERP\nLoad"),] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 3) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_manual(values = my.cols) + ylim(0.0, NA) 
P$layers <- c(geom_boxplot(outlier.shape = NA), P$layers)
P <- P + theme(legend.position = "bottom", legend.key.size = unit(1, 'cm'), axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank()) 
P <- P + guides(fill=guide_legend(nrow=1,byrow=TRUE, label.position = "bottom"))
my.legend <- get_legend(P)
P <- P + theme(legend.position = "none")


Q <- output[output$Score == "BAD_Mutations\nLoad",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 3) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_manual(values = my.cols)  + 
  scale_y_continuous(position = "right", limits = c(0, NA)) 
Q$layers <- c(geom_boxplot(outlier.shape = NA), Q$layers)
Q <- Q + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank()) 

R <- output[c(output$Score == "Synonymous\nMutations" | output$Score == "Nonsynonymous\nMutations"),] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 3) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_manual(values = my.cols)  + ylim(0.0, NA) 
R$layers <- c(geom_boxplot(outlier.shape = NA), R$layers)
R <- R + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank()) 


#S <- output[output$Score == "Nonsynonymous\nMutations",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 1.5) +
  #facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_brewer(palette = "Accent") + ylim(0.0, NA) 
#S$layers <- c(geom_boxplot(outlier.shape = NA), S$layers)
#S <- S + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank())  





#U <- output[output$Score == "All Genes",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 1.5) +
#  facet_grid(Phylogeny ~ Score, scales = "free") + ylab("% BAD_Mutations / Nonsynonymous Mutations") +
#  theme_bw() + scale_fill_brewer(palette = "Accent") #+ ylim(0.0, NA) 
#U$layers <- c(geom_boxplot(outlier.shape = NA), U$layers)
#U <- U + theme(legend.position = "none", axis.title.x = element_blank(), strip.text.y = element_blank())  

#####Homoeologs 

homoeolog_loads <- read.table("Homoeologs_Load_per_Pop.NO_GC.txt", header = TRUE)
new_sums <- data.frame(Sample = names(homoeolog_loads[,c(3:258)]), Load = colSums(homoeolog_loads[,c(3:258)]))
new_sums$Species  = rep(c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))
new_sums$Phylogeny = rep(c("Basal", "Diploid", "Polyploid", "Seg"), each = 64)
new_sums$Score = rep(c(rep("GERP\nLoad", times = 16), rep("BAD_Mutations\nLoad", times = 16), rep("Synonymous\nMutations", times=16), rep("Nonsynonymous\nMutations", times=16)))
new_sums$Subgenome = rep(c("At", "Dt"), each = 8)
new_sums$Scale = rep("Homoeologs\n(no GC)")
new_sums <- new_sums[!new_sums$Phylogeny == "Seg",]

#output <- rbind(new_sums, BAD)
new_sums <- new_sums[!c(new_sums$Species == "A1" & new_sums$Subgenome == "Dt" & new_sums$Phylogeny != "Basal"),]
new_sums <- new_sums[!c(new_sums$Species == "D5" & new_sums$Subgenome == "At" & new_sums$Phylogeny != "Basal"),]
new_sums <- new_sums[!c(new_sums$Species == "A1" & new_sums$Phylogeny == "Polyploid"),]
new_sums <- new_sums[!c(new_sums$Species == "D5" & new_sums$Phylogeny == "Polyploid"),]
#output <- output[!c(output$Phylogeny == "Segregating"),]

new_sums$Score <- factor(new_sums$Score, levels = c("Synonymous\nMutations", "Nonsynonymous\nMutations", "GERP\nLoad", "BAD_Mutations\nLoad"))


BAD <- new_sums[new_sums$Score == "BAD_Mutations\nLoad",]
Nonsynon <- new_sums[new_sums$Score == "Nonsynonymous\nMutations",]
BAD$percent <- (BAD$Load / Nonsynon$Load) * 100
BAD$Load <- BAD$percent
BAD <- BAD[,-c(1,8)]
BAD$Score <- rep("Homoeologs")

Percent_Bad <- rbind(Percent_Bad, BAD)
#Percent_Bad <- BAD

P2 <- new_sums[c(new_sums$Score == "GERP\nLoad"),] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 3) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_manual(values = my.cols) + ylim(0.0, NA) 
P2$layers <- c(geom_boxplot(outlier.shape = NA), P2$layers)
P2 <- P2 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank()) 

Q2 <- new_sums[new_sums$Score == "BAD_Mutations\nLoad",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 3) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_manual(values = my.cols) + 
  scale_y_continuous(position = "right", limits = c(0, NA)) 
Q2$layers <- c(geom_boxplot(outlier.shape = NA), Q2$layers)
Q2 <- Q2 + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank())
Q2 <- Q2 + theme(legend.position = "none")


R2 <- new_sums[c(new_sums$Score == "Synonymous\nMutations" | new_sums$Score == "Nonsynonymous\nMutations"),] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 3) +
  facet_grid(Phylogeny ~ Score, scales = "free") + theme_bw() + scale_fill_manual(values = my.cols)  + ylim(0.0, NA)  
R2$layers <- c(geom_boxplot(outlier.shape = NA), R2$layers)
R2 <- R2 + theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), strip.text.y = element_blank()) 



#png("Figure_Loads_Mutations_per_Phylogeny.png", pointsize = 12, width = 4500, height = 3000, res = 600)
(R2 + P2 + Q2 + plot_layout(widths = c(2.3,1,1.1)))/ my.legend  + plot_layout(heights = c(2,0.25))
#dev.off()
ggsave("Figure_Loads_Mutations_per_Phylogeny.pdf", height = 12, width = 15, units = "cm")


#png("Figure_Loads_Mutations_per_Phylogeny_full.png", pointsize = 12, width = 4500, height = 5500, res = 600)
(R + P + Q + plot_layout(widths = c(2.3,1,1.1)))  / my.legend / (R2 + P2 + Q2 + plot_layout(widths = c(2.3,1,1.1)))  + plot_layout(heights = c(2,0.25,2))
#dev.off()
ggsave("Figure_Loads_Mutations_per_Phylogeny_full.pdf", height = 22, width = 15, units = "cm")


X <- Percent_Bad %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 2) +
  facet_grid(Phylogeny ~ Score, scales = "free") + ylab("% BAD_Mutations / Nonsynonymous Mutations") +
  theme_bw() + scale_fill_manual(values = my.cols) #+ ylim(0.0, NA) 
X$layers <- c(geom_boxplot(outlier.shape = NA), X$layers)
X <- X + theme(axis.title.x = element_blank(), strip.text.y = element_blank())  

#png("Figure_Loads_Mutations_percentages.png", pointsize = 12, width = 3000, height = 3000, res = 600)
X
#dev.off()
ggsave("Figure_Loads_Mutations_percentages_full.pdf", height = 11, width = 8, units = "cm")

X <- Percent_Bad[Percent_Bad$Score == "Homoeologs",] %>% ggdotplot(x="Subgenome", y="Load", fill = "Species", size = 2) +
  facet_grid(Phylogeny ~ Score, scales = "free") + ylab("% BAD_Mutations / Nonsynonymous Mutations") +
  theme_bw() + scale_fill_manual(values = my.cols) #+ ylim(0.0, NA) 
X$layers <- c(geom_boxplot(outlier.shape = NA), X$layers)
X <- X + theme(axis.title.x = element_blank(), strip.text.y = element_blank()) 
X
ggsave("Figure_Loads_Mutations_percentages.pdf", height = 11, width = 8, units = "cm")

####Not Used From Here Down 
homoeolog_loads <- read.table("../Homoeologs_Load_per_Pop.txt", header = TRUE)
bed <- read.table("../pSONIC.homoeologs.CDS.correctCHR.correctENDPOS.orientation.bed")
names(bed) <- c("Chr", "Start", "End", "Orientation", "Gene")
bed$Gene <- substr(bed$Gene,1,nchar(bed$Gene)-2)
bed$Gene <- gsub('\\.', '', bed$Gene)

find_length <- function(gene){
  temp <- bed[bed$Gene == gene,]
  return(sum(temp$End - temp$Start))
}

find_beginning <- function(gene){
  temp <- bed[bed$Gene == gene,]
  return(min(temp$Start))
}

find_chrom <- function(gene){
  temp <- bed[bed$Gene == gene,]
  return(min(temp$Chr))
}


homoeolog_loads$At_Length <- unlist(lapply(homoeolog_loads$At, find_length))
homoeolog_loads$Dt_Length <- unlist(lapply(homoeolog_loads$Dt, find_length))
homoeolog_loads$At_POS <- unlist(lapply(homoeolog_loads$At, find_beginning))
homoeolog_loads$Dt_POS <- unlist(lapply(homoeolog_loads$Dt, find_beginning))
homoeolog_loads$At_CHROM <- unlist(lapply(homoeolog_loads$At, find_chrom))
homoeolog_loads$Dt_CHROM <- unlist(lapply(homoeolog_loads$Dt, find_chrom))


#homoeolog_pos <- homoeolog_loads[,c(1,2,)]
#homoeolog_pos$Name = substr(homoeolog_pos$Name,1,nchar(homoeolog_pos$Name)-2)
#homoeolog_pos$Name = gsub('\\.', '', homoeolog_pos$Name)

#pairs.genes <- homoeolog_loads[,c(1,2)]

#pairs.genes <- merge(pairs.genes, homoeolog_pos, by.x="At", by.y="Name", all.x = TRUE)
#names(pairs.genes) <- c("At", "Dt", "At_CHROM", "At_POS")
#pairs.genes <- merge(pairs.genes, homoeolog_pos, by.x="Dt", by.y="Name", all.x = TRUE)
#names(pairs.genes) <- c("Dt", "At", "At_CHROM", "At_POS", "Dt_CHROM", "Dt_POS")
homoeolog_loads <- homoeolog_loads[order(homoeolog_loads$At_POS),]
homoeolog_loads$At_order <- ave(homoeolog_loads$At_POS, homoeolog_loads$At_CHROM, FUN=seq_along)
homoeolog_loads <- homoeolog_loads[order(homoeolog_loads$Dt_POS),]
homoeolog_loads$Dt_order <- ave(homoeolog_loads$Dt_POS, homoeolog_loads$Dt_CHROM, FUN=seq_along)

#all.scores <- merge(pairs.genes, homoeolog_loads, by.x = c("At", "Dt"), by.y = c("At", "Dt"))
#GERP <- homoeolog_loads[,c("At", "Dt", "At_CHROM", "At_order", "Dt_CHROM", "At_Length", "Dt_Length",
#                     "At_AD1_BAD_load_dip", "At_AD3_BAD_load_dip", "At_AD4_BAD_load_dip", "At_AD5_BAD_load_dip", "At_AD6_BAD_load_dip", "At_AD7_BAD_load_dip", "At_A1_BAD_load_dip", "At_D5_BAD_load_dip",
#                     "Dt_AD1_BAD_load_dip", "Dt_AD3_BAD_load_dip", "Dt_AD4_BAD_load_dip", "Dt_AD5_BAD_load_dip", "Dt_AD6_BAD_load_dip", "Dt_AD7_BAD_load_dip", "Dt_A1_BAD_load_dip", "Dt_D5_BAD_load_dip")]
GERP <- homoeolog_loads[,c("At", "Dt", "At_CHROM", "At_order", "Dt_CHROM", "At_Length", "Dt_Length",
                     "At_AD1_BAD_load", "At_AD3_BAD_load", "At_AD4_BAD_load", "At_AD5_BAD_load", "At_AD6_BAD_load", "At_AD7_BAD_load", "At_A1_BAD_load", "At_D5_BAD_load",
                     "Dt_AD1_BAD_load", "Dt_AD3_BAD_load", "Dt_AD4_BAD_load", "Dt_AD5_BAD_load", "Dt_AD6_BAD_load", "Dt_AD7_BAD_load", "Dt_A1_BAD_load", "Dt_D5_BAD_load")]

GERP$Score <- rep("GERP")
#GERPA <- GERP[,c("At", "At_CHROM", "At_Length", "At_order","At_AD1_BAD_load_dip", "At_AD3_BAD_load_dip", "At_AD4_BAD_load_dip", "At_AD5_BAD_load_dip", "At_AD6_BAD_load_dip", "At_AD7_BAD_load_dip", "At_A1_BAD_load_dip", "At_D5_BAD_load_dip", "Score")]
#GERPD <- GERP[,c("Dt", "At_CHROM", "Dt_Length", "At_order","Dt_AD1_BAD_load_dip", "Dt_AD3_BAD_load_dip", "Dt_AD4_BAD_load_dip", "Dt_AD5_BAD_load_dip", "Dt_AD6_BAD_load_dip", "Dt_AD7_BAD_load_dip", "Dt_A1_BAD_load_dip", "Dt_D5_BAD_load_dip", "Score")]
GERPA <- GERP[,c("At", "At_CHROM", "At_Length", "At_order","At_AD1_BAD_load", "At_AD3_BAD_load", "At_AD4_BAD_load", "At_AD5_BAD_load", "At_AD6_BAD_load", "At_AD7_BAD_load", "At_A1_BAD_load", "At_D5_BAD_load", "Score")]
GERPD <- GERP[,c("Dt", "At_CHROM", "Dt_Length", "At_order","Dt_AD1_BAD_load", "Dt_AD3_BAD_load", "Dt_AD4_BAD_load", "Dt_AD5_BAD_load", "Dt_AD6_BAD_load", "Dt_AD7_BAD_load", "Dt_A1_BAD_load", "Dt_D5_BAD_load", "Score")]

names(GERPA) <- c("Gene", "CHROM", "Length", "order", "AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5", "Score")
names(GERPD) <- c("Gene", "CHROM", "Length", "order", "AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5", "Score")
GERPA <- GERPA[with(GERPA, order(CHROM, order)),]
GERPD <- GERPD[with(GERPD, order(CHROM, order)),]
GERPD$CHROM <- GERPD$CHROM + 13
polyploids <- c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7")
GERP <- rbind(GERPA, GERPD)

find_roll_ave <- function(CHR, Species, number, step){
  temp <- GERP[GERP$CHROM == CHR,]
  temp <- temp[order(temp$order),]
  spec1 <- Species[1]
  result1 <- rollapply(temp[[spec1]], width=number, by=step, FUN=mean, align="left")
  result <- data.frame(CHROM = rep(CHR, length(result1)), order = seq(1:length(result1)))
  result[, spec1] <- result1
  Species <- Species[-1]
  for(i in Species){
    result[[i]] <- rollapply(temp[[i]], width=number, by=step, FUN=mean, align="left")
  }
  #sapply(Species, rollapply, width=5, by=2, FUN=mean, align="left")
  return(result)
}

output <- find_roll_ave(1, c('AD1', 'AD3', 'AD4', 'AD5', 'AD6', 'AD7', 'A1', 'D5', 'Length'), 25, 1)
for(i in c(2:26)){
  temp.out <- find_roll_ave(i, c('AD1', 'AD3', 'AD4', 'AD5', 'AD6', 'AD7', 'A1', 'D5', 'Length'), 25, 1)
  output <- rbind(output, temp.out)
}


GERPA <- output[output$CHROM < 14,]
GERPD <- output[output$CHROM > 13,]
GERPA$Subgenome <- rep("At")
GERPD$Subgenome <- rep("Dt")
GERPD$CHROM <- GERPD$CHROM - 13
GERPA$D5 <- rep(0)
GERPD$A1 <- rep(0)


GERPTotal <- data.frame(CHROM = GERPA$CHROM, order = GERPA$order, Length = GERPA$Length, Subgenome = rep("Total"))

for(i in polyploids){
  GERPA[[i]] <- GERPA[[i]] / GERPA$D5
  GERPD[[i]] <- GERPD[[i]] / GERPD$A1
}

GERPA$A1 <- GERPA$A1 / GERPA$Length
GERPD$D5 <- GERPD$D5 / GERPD$Length

for(i in polyploids){
  GERPA[[i]] <- GERPA[[i]] - GERPA$A1
  GERPD[[i]] <- GERPD[[i]] - GERPD$D5
}

for(i in polyploids){
  GERPA[[i]] <- GERPA[[i]] / GERPA$Length
  GERPD[[i]] <- GERPD[[i]] / GERPD$Length
}


for(i in polyploids){
  GERPA[[i]] <- log2(GERPA[[i]])
  GERPD[[i]] <- log2(GERPD[[i]])
}

for(i in polyploids){
  GERPTotal[[i]] <- (GERPA[[i]] + GERPD[[i]])
}

GERPTotal$A1 <- GERPA$A1 + GERPD$A1
GERPTotal$D5 <- GERPA$D5 + GERPD$D5
GERP <- rbind(GERPA, GERPD)#, GERPTotal)
GERP$Diploid <- GERP$A1 + GERP$D5



GERP_melt <- melt(GERP, measure.vars = c("Diploid"), #c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7"),
                    variable.name = "Species", value.name = "Load")
#GERP_melt <- GERP_melt[complete.cases(GERP_melt),]

#ave(GERP_melt$Load, GERP_melt$Species, FUN = mean)

my.cols <- my.cols[2:7]


png("Figure_Loads_Chromosome_basal_only_diploids_25_1.png", pointsize = 12, width = 4000, height = 1000, res = 600)
ggplot(GERP_melt[c(GERP_melt$CHROM == 3),], aes(x=order, y=Load, color = Subgenome)) + 
  geom_point(size = 0.01) + theme_bw() + guides(colour = guide_legend(override.aes = list(size=10))) + 
  facet_grid(Species~., scales = "free") + scale_color_manual(values = my.cols)
dev.off()

hist(GERPA$AD1 - GERPD$AD1, breaks = 100)

summary(GERPA$AD1 - GERPD$AD1)

sd(GERPA$AD1)
sd(GERPD$AD1)
sd(GERPTotal$AD1)

ggplot(pairs.genes, aes(x=Dt_POS, y=Dt_order)) + 
  geom_point() + facet_wrap(.~Dt_CHROM, scales = "free")



for(i in polyploids){
  GERPA[[i]] <- GERPA[[i]] / GERPA$A1
  GERPD[[i]] <- GERPD[[i]] / GERPD$D5
}

for(i in polyploids){
  GERPA[[i]] <- GERPA[[i]] / GERPA$Length
  GERPD[[i]] <- GERPD[[i]] / GERPD$Length
}
for(i in polyploids){
  GERPA[[i]] <- log2(GERPA[[i]])
  GERPD[[i]] <- log2(GERPD[[i]])
}

plot(GERPA$AD1 , GERPD$AD1)
abline(lm(GERPA$AD1 ~ GERPD$AD1))
