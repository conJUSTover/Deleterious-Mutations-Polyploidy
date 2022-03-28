library(stringr)
library(ggplot2)
library(ggupset)
library(tidyverse)
library(scales)

setwd("~/Documents/GERP/pSONIC_files/map_D_no_HE/sort_SNPs/")

snps <- read.table("GERP06_BAD.frqs.pSONIC.NO_GC.txt", header = TRUE)

classify_frq <- function(SNP){
  if(is.na(SNP)){
    return(NA)
  } else if(SNP == 0){
    return("W")
  } else if(SNP == 1){
    return("B")
  } else {
    return("G")
  }
}

snps$A1U <- sapply(snps$A1, classify_frq)
snps$D5U <- sapply(snps$D5, classify_frq)
snps$AD1U <- sapply(snps$AD1, classify_frq)
snps$AD3U <- sapply(snps$AD3, classify_frq)
snps$AD4U <- sapply(snps$AD4, classify_frq)
snps$AD5U <- sapply(snps$AD5, classify_frq)
snps$AD6U <- sapply(snps$AD6, classify_frq)
snps$AD7U <- sapply(snps$AD7, classify_frq)

At <- snps[snps$CHROM < 14,]
Dt <- snps[snps$CHROM > 13,]

At <- At[complete.cases(At[,c("A1U", "AD4U", "AD5U", "AD3U", "AD1U", "AD6U", "AD7U")]),]
Dt <- At[complete.cases(At[,c("D5U", "AD4U", "AD5U", "AD3U", "AD1U", "AD6U", "AD7U")]),]

At$UpSet <- paste(At$A1U, At$AD4U, At$AD5U, At$AD3U, At$AD1U, At$AD6U, At$AD7U, sep = "")
Dt$UpSet <- paste(Dt$D5U, Dt$AD4U, Dt$AD5U, Dt$AD3U, Dt$AD1U, Dt$AD6U, Dt$AD7U, sep = "")

At %>%
  count(UpSet) -> At2

Dt %>%
  count(UpSet) -> Dt2

#At2 <- At2[At2$n > 100,]

#ggplot(At2, aes(x = reorder(UpSet, -n), y = n)) + 
#  geom_bar(stat = "identity") + 
#  theme_bw() + 
# theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1), 
#        axis.title.x = element_blank()) + 
#  ylab("Number of Homoeologous SNP Loci")


find_combinations <- function(W, L){
  x <- expand.grid(rep(list(c('W', 'B', 'G')), L-W))
  y <- do.call(paste0, x)
  bad_strings <- c()
  for(i in y){
    B <- lengths(regmatches(i, gregexpr("B", i)))
    G <- lengths(regmatches(i, gregexpr("G", i)))
    if(B == 1 & G == 0){
      bad_strings <- c(bad_strings, i)
    } else if(B == 0 & G == 0){
      bad_strings <- c(bad_strings, i)
    } else if(B == 0 & G == 1){
      bad_strings <- c(bad_strings, i)
    } else if(substr(i,1,1) == "W"){
      bad_strings <- c(bad_strings, i)
    } else if(G == 0 & B > 0){
      result <- substr(i,(nchar(i)+1)-B,nchar(i))
      if(result == paste(rep("B", B), collapse="")){
        bad_strings <- c(bad_strings, i)
      }
    } 
  }
  y <- y[!y %in% bad_strings]
  y <- paste(paste0(rep("W", W), collapse=""), y, sep="")
  return(y)
}


phylogenetic_nonhomoplasies <- c("BBBBBBB", "WBBBBBB","WWBBBBB","WWWBBBB","WWWWBBB","WWWWWBB","WWWWWWB","WWWWWBW","WWWWBWW","WWWBWWW","WWBWWWW","WBWWWWW","BWWWWWW")
terminal_segs <- c("WWWWWWG", "WWWWWGW","WWWWGWW","WWWGWWW","WWGWWWW","WGWWWWW","GWWWWWW")
AD67_phylo_homos <- find_combinations(5,7)
AD167_phylo_homos<- find_combinations(4,7)
AD3167_phylo_homos<- find_combinations(3,7)
AD53167_phylo_homos<- find_combinations(2,7)
AD453167_phylo_homos<- find_combinations(1,7)
ADA453167_phylo_homos<- find_combinations(0,7)


AT_NoHomoplasies <- At2[At2$UpSet %in% phylogenetic_nonhomoplasies,]
At_terminal_segs <- At2[At2$UpSet %in% terminal_segs,]
A67_homos <- sum(At2[At2$UpSet %in% AD67_phylo_homos,"n"])
A167_homos <- sum(At2[At2$UpSet %in% AD167_phylo_homos,"n"])
A3167_homos <- sum(At2[At2$UpSet %in% AD3167_phylo_homos,"n"])
A53167_homos <- sum(At2[At2$UpSet %in% AD53167_phylo_homos,"n"])
A453167_homos <- sum(At2[At2$UpSet %in% AD453167_phylo_homos,"n"])
AA453167_homos <- sum(At2[At2$UpSet %in% ADA453167_phylo_homos,"n"])

DT_NoHomoplasies <- Dt2[Dt2$UpSet %in% phylogenetic_nonhomoplasies,]
Dt_terminal_segs <- Dt2[Dt2$UpSet %in% terminal_segs,]
D67_homos <- sum(Dt2[Dt2$UpSet %in% AD67_phylo_homos,"n"])
D167_homos <- sum(Dt2[Dt2$UpSet %in% AD167_phylo_homos,"n"])
D3167_homos <- sum(Dt2[Dt2$UpSet %in% AD3167_phylo_homos,"n"])
D53167_homos <- sum(Dt2[Dt2$UpSet %in% AD53167_phylo_homos,"n"])
D453167_homos <- sum(Dt2[Dt2$UpSet %in% AD453167_phylo_homos,"n"])
DD453167_homos <- sum(Dt2[Dt2$UpSet %in% ADA453167_phylo_homos,"n"])


