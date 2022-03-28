library(reshape2)
library(RColorBrewer)
library(ggplot2)
setwd("~/Documents/GERP/good_filtering_results/")

infile <- read.table("GERP06_BAD.frqs.txt", header = TRUE)

my.cols <- brewer.pal(9, "Set1")
my.cols[9] <- "#000000"
my.cols <- my.cols[c(1,2,4:9)]

bonferonni <- 0.05 / 967154
  
GERP_bypos <- function(df, sp){
  return(sum(df[,sp] * df$GERP, na.rm = TRUE))
}

GERP_bypos_all <- function(start, window, CHR, Species){
  temp <- infile[c(infile$GERP > 0 & !is.na(infile$GERP) & infile$POS < start & infile$POS > (start - window) & infile$CHROM == CHR),]
  output <- sapply(Species, GERP_bypos, df=temp)
  return(output)
}

BAD_bypos <- function(df, sp){
  return(sum(df[,sp], na.rm = TRUE))
}

BAD_bypos_all <- function(start, window, CHR, Species){
  temp <- infile[c(infile$MaskedP.value < bonferonni & !is.na(infile$MaskedP.value) & infile$POS < start & infile$POS > (start - window) & infile$CHROM == CHR),]
  output <- sapply(Species, BAD_bypos, df=temp)
  return(output)
}

At <- infile[infile$CHROM < 14,]
Dt <- infile[infile$CHROM > 13,]
At$average <- rowMeans(At[,c(5:10,12)], na.rm = TRUE)
Dt$average <- rowMeans(Dt[,c(5:11)], na.rm = TRUE)
At <- At[c(At$average > 0 & At$average < 1 & !is.na(At$average)),]
Dt <- Dt[c(Dt$average > 0 & Dt$average < 1 & !is.na(Dt$average)),]
infile <- rbind(At, Dt)

windows <- 10000000 #100 KB
steps <- 300000 #30KB 

find_load_segment <- function(windows, steps, chrom){
 ch01 <- infile[infile$CHROM == chrom,]
 max.pos <- max(ch01$POS)
 max.round <- ceiling(max.pos/steps) * steps
 output_BAD <- data.frame(CHROM = rep(chrom), Steps = seq(windows, max.round, by=steps))
 output_GERP <- data.frame(CHROM = rep(chrom), Steps = seq(windows, max.round, by=steps))
 
 output_BAD[, c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5")] <- t(apply(output_BAD, 1, function(x) BAD_bypos_all(x[2], windows, chrom, c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
 output_GERP[, c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5")] <- t(apply(output_GERP, 1, function(x) GERP_bypos_all(x[2], windows, chrom, c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "A1", "D5"))))
 output_BAD$Score <- rep("BAD_Mutations")
 output_GERP$Score <- rep("GERP_Load")
  output <- rbind(output_BAD, output_GERP)
  return(output)
}

results <- find_load_segment(windows, steps, 1)

for(i in c(2:26)){
  results <- rbind(results, find_load_segment(windows, steps, i))
  print(i)
}

At <- results[results$CHROM < 14,]
At$Subgenome <- rep("At")

Dt <- results[results$CHROM > 13,]
Dt$Subgenome <- rep("Dt")
Dt$CHROM <- Dt$CHROM - 13

results <- rbind(At, Dt)

final <- melt(results, measure.vars = c("A1", "AD1", "AD3", "AD4", "AD5", "AD6", "AD7", "D5"),
              variable.name = "Species", value.name = "Load")
final <- final[!c(final$Subgenome == "At" & final$Species == "D5"),]
final <- final[!c(final$Subgenome == "Dt" & final$Species == "A1"),]


png("Figure_Loads_Chromosome.png", pointsize = 12, width = 4000, height = 7000, res = 600)
ggplot(final[final$Score == "BAD_Mutations",], aes(x=Steps, y=Load, color = Species)) + 
  geom_point(size = 0.01) + theme_bw() + guides(colour = guide_legend(override.aes = list(size=10))) + 
  facet_grid(CHROM~Subgenome, scales = "free") + scale_color_manual(values = my.cols)
dev.off()

At$AD1 <- At$AD1 / At$A1
At$AD3 <- At$AD3 / At$A1
At$AD4 <- At$AD4 / At$A1
At$AD5 <- At$AD5 / At$A1
At$AD6 <- At$AD6 / At$A1
At$AD7 <- At$AD7 / At$A1

Dt$AD1 <- Dt$AD1 / Dt$D5
Dt$AD3 <- Dt$AD3 / Dt$D5
Dt$AD4 <- Dt$AD4 / Dt$D5
Dt$AD5 <- Dt$AD5 / Dt$D5
Dt$AD6 <- Dt$AD6 / Dt$D5
Dt$AD7 <- Dt$AD7 / Dt$D5

results <- rbind(At, Dt)
final_again <- melt(results, measure.vars = c("AD1", "AD3", "AD4", "AD5", "AD6", "AD7"),
                          variable.name = "Species", value.name = "Load")
#final_again <- final_again[c(final_again$Species != "A1" & final_again$Species != "D5"),]

my.cols <- my.cols[2:7]
png("Figure_Loads_Chromosome_relative_increase.png", pointsize = 12, width = 4000, height = 7000, res = 600)
ggplot(final_again[final_again$Score == "BAD_Mutations",], aes(x=Steps, y=Load, color = Species)) + 
  geom_point(size = 0.01) + theme_bw() + guides(colour = guide_legend(override.aes = list(size=10))) + 
  facet_grid(CHROM~Subgenome, scales = "free") + scale_color_manual(values = my.cols)
dev.off()

