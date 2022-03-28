
setwd("/work/LAS/jfw-lab/del_mut_gvcf/redo_f2")

read_counts_file <- function(sample, outfile){
  infile <- read.table(paste(sample, ".sort.coverage.txt", sep=''))
  print(sample)
  names(infile) <- c("Chr", "Start", "Stop", "GeneID", "GenePos", "Depth")
  output <- aggregate(infile$Depth, by=list(infile$Chr, infile$Start, infile$Stop, infile$GeneID), FUN=sum)
  rm(infile)
  names(output) <- c("Chr", "Start", "Stop", "GeneID", "ReadDepth")
  output$Sample <- rep(sample)
  output <- rbind(output, outfile)
  return(output)
}

infile <- read.table("100.sort.coverage.txt")
print("100 read")
names(infile) <- c("Chr", "Start", "Stop", "GeneID", "GenePos", "Depth")

output <- aggregate(infile$Depth, by=list(infile$Chr, infile$Start, infile$Stop, infile$GeneID), FUN=sum)
rm(infile)
names(output) <- c("Chr", "Start", "Stop", "GeneID", "ReadDepth")
output$Sample <- rep("100")

other_samples <- c("101","102","103","104","105","106","107","108","109","301","302","303","304","305","306","307","401","402","403","405","406","502","504","505","507","601","602","603","604","605","606","701","702")

for (i in other_samples){
  output <- read_counts_file(i, output)
}

write.table(output, "F260.ReadCounts.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

