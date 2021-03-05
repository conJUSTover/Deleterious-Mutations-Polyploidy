setwd("/work/LAS/jfw-lab/jconover/deleterious_mutations/BAD_runs/redo_Australian")

BAD_table <- read.table("Combined_Report.txt", header = TRUE)

subs <- read.table("all.GT.subs.longform.txt")
names(subs) <- c("gene", "aapos", "altaa", "SNP Name")

##Create new variable for gene_aaPOS
BAD_table$GeneAA <- apply(BAD_table[,c("GeneID", "CDSPos")], 1, paste, sep='', collapse = "_")
subs$GeneAA <- apply(subs[,c("gene", "aapos")], 1, paste, sep='', collapse = "_")

##Merge SNP reference and BAD_table, keeping all sub entries and duplicating 
##SNPs occuring in same codon 
new_BAD <- merge(BAD_table, subs, by.x = "GeneAA", by.y = "GeneAA", all.y = TRUE)
new_BAD$GeneAA <- gsub(" ",  "", new_BAD$GeneAA)

#write to new file
write.table(new_BAD, file = "all.GT.merged_results.txt", quote = FALSE, sep = "\t", row.names = FALSE)

