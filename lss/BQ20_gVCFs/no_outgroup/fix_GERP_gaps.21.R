setwd("/work/LAS/jfw-lab/del_mut_gvcf/BQ20_gVCFs/vcfs")


fasta <- read.table('roast.21.msa.Ghirsutum.seqlines', header=F)
fasta$rows <- as.numeric(rownames(fasta))
rates <- read.table('roast.21.msa.in.rates', header=F)
rates$rowname <- as.numeric(rownames(rates))

indices_to_add <- data.frame(x=fasta$rows[fasta$V1 == 'N'])
indices_to_add$row <- as.numeric(rownames(indices_to_add))
indices_to_add$newindex <- indices_to_add$x + 0.5 - indices_to_add$row

add_rows = data.frame(V1 = c(rep(0, length(indices_to_add$newindex))), 
                      V2 = c(rep(0, length(indices_to_add$newindex))),
                      rowname = indices_to_add$newindex)

rates <- rbind(rates, add_rows)
rates <- rates[order(rates$rowname),]
rates <- rates[c("V1", "V2")]
rownames(rates) <- NULL

write.table(rates, "roast.21.msa.in.new.rates", sep='\t', row.names=F, col.names = F)

rates$rowname <- as.numeric(rownames(rates))


#Find gaps in all the other seqs 
fasta$Gk <- read.table('roast.21.msa.Gkirkii.seqlines', header=F)
fasta$Tc <- read.table('roast.21.msa.Tcacao.seqlines', header=F)
fasta$Cc <- read.table('roast.21.msa.Cclementina.seqlines', header=F)
fasta$Cm <- read.table('roast.21.msa.Cmaxima.seqlines', header=F)
fasta$Cp <- read.table('roast.21.msa.Cpapaya.seqlines', header=F)
fasta$Cs <- read.table('roast.21.msa.Csinensis.seqlines', header=F)
fasta$Fv <- read.table('roast.21.msa.Fvesca.seqlines', header=F)
fasta$Pp <- read.table('roast.21.msa.Ppersica.seqlines', header=F)
fasta$Vv <- read.table('roast.21.msa.Vvinifera.seqlines', header=F)

gaps <- as.numeric(rownames(fasta[fasta$Gk == "-" | fasta$Tc == "-" | 
                                    fasta$Cc == "-" | fasta$Cm == "-" | 
                                    fasta$Cp == "-" | fasta$Cs == "-" | 
                                    fasta$Fv == "-" | fasta$Pp == "-" | 
                                    fasta$Vv == "-",]))

rates <- rates[! rates$rowname %in% gaps,]

gap_readd = data.frame(V1 = c(rep(0, length(gaps))), 
                      V2 = c(rep(0, length(gaps))),
                      rowname = gaps)

rates <- rbind(rates, gap_readd)
rates <- rates[order(rates$rowname),]
rates <- rates[c("V1", "V2")]
rownames(rates) <- NULL

write.table(rates, "roast.21.msa.in.gaps.rates", sep='\t', row.names=F, col.names = F)
