vcf <- read.table(gzfile('all.frq.txt.gz'), header = T)
print("File loaded successfully")

REF <- vcf[vcf$X090 == "0/0",]
ALT <- vcf[vcf$X090 == "1/1",]

print("vcf Split into REF and ALT successfully")

names(REF) <- c("Chr","POS","090","Anc_101","Der_101","Filter_101","Anc_10","Der_10","Filter_10","Anc_11","Der_11","Filter_11","Anc_12","Der_12","Filter_12","Anc_13","Der_13","Filter_13","Anc_21","Der_21","Filter_21","Anc_22","Der_22","Filter_22", "Anc_23","Der_23","Filter_23","Anc_30","Der_30","Filter_30","Anc_40","Der_40","Filter_40","Anc_50","Der_50","Filter_50","Anc_60","Der_60","Filter_60","Anc_70","Der_70","Filter_70","Anc_80","Der_80","Filter_80","Anc_90","Der_90","Filter_90","Anc_91","Der_91","Filter_91","Gen_REF","Gen_Anc","GERPgap0","GERP")
                
print("Name of REF changed successfully")

names(ALT) <- c("Chr","POS","090","Der_101","Anc_101","Filter_101","Der_10","Anc_10","Filter_10","Der_11","Anc_11","Filter_11","Der_12","Anc_12","Filter_12","Der_13","Anc_13","Filter_13","Der_21","Anc_21","Filter_21","Der_22","Anc_22","Filter_22", "Der_23","Anc_23","Filter_23","Der_30","Anc_30","Filter_30","Der_40","Anc_40","Filter_40","Der_50","Anc_50","Filter_50","Der_60","Anc_60","Filter_60","Der_70","Anc_70","Filter_70","Der_80","Anc_80","Filter_80","Der_90","Anc_90","Filter_90","Der_91","Anc_91","Filter_91","Gen_REF","Gen_Anc","GERPgap0","GERP")

print("Nanes of ALT changed successfully")

vcf <- rbind(REF, ALT)

print("Rbind success")

vcf <- vcf[,c("Chr","POS","090","Der_101","Filter_101","Der_10","Filter_10","Der_11","Filter_11","Der_12","Filter_12","Der_13","Filter_13","Der_21","Filter_21","Der_22","Filter_22", "Der_23","Filter_23","Der_30","Filter_30","Der_40","Filter_40","Der_50","Filter_50","Der_60","Filter_60","Der_70","Filter_70","Der_80","Filter_80","Der_90","Filter_90","Der_91","Filter_91","Gen_REF","GERP")]

print("vcf pruned to only include DAF")

#mask <- vcf[,seq(5, ncol(vcf) -2, 2)] == "PASS"

#print("MASK created")

#vcf[,seq(5, ncol(vcf) -2, 2)] <- mask

#print("MASK applied")

write.table(vcf, "all.frq.DER.txt", sep='\t', row.names=F, col.names = T, quote=F)

print("first table written")

vcf <- vcf[vcf$GERP > 0,]

print("Table Pruning GERP > 0 done")

write.table(vcf, "all.frq.DER.GERP06.txt", sep='\t', row.names=F, col.names = T, quote=F)
