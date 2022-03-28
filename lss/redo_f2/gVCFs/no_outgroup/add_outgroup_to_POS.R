library(stringr)
library(dplyr)
library(reshape2)

infile <- read.table("pSONIC.homoeoSNPs.GT.pos", header = FALSE)
names(infile) <- c("CHROM", "POS")

infile$CHROM <- str_pad(infile$CHROM, 2, pad = "0")
infile$POS <- as.character(infile$POS)
infile$CHROM_POS <-  apply(infile[,c("CHROM", "POS")], 1, paste, collapse = "_")
infile$POS <- as.integer(infile$POS)

GERP_vcf <- infile[FALSE,]

for(i in c(01:26)){
  j <- str_pad(i, 2, pad = "0")
  print(j)
  GERP_name <- paste("../GERP_SCORES/roast", j, "msa.in.Tcacao.seqlines", sep = ".")
  GERP <- read.table(GERP_name, header=FALSE)
#  names(GERP) <- c('Neutral Rate','RS Score')
  temp_frqs <- infile[infile$CHROM == j,]
  temp_frqs$ANC <- GERP[temp_frqs$POS,1]
  GERP_vcf <- dplyr::bind_rows(temp_frqs, GERP_vcf)
}

write.table(GERP_vcf, file = "pSONIC.homoeoSNPs.GT.pos.Tcacao.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
