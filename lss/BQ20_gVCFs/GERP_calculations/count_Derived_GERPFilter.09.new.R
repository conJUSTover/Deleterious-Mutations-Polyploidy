library(data.table)
setwd("/work/LAS/jfw-lab/del_mut_gvcf/BQ20_gVCFs/GERP_calculations")

vcf_header <- c('#09','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','090','1010','1011','1012','1013','1014','1015','1016','1017','1018','1019','100','101','102','103','104','105','106','107','108','109','110','111','112','113','114','115','116','117','118','119','120','121','122','123','124','125','126','127','128','129','130','131','132','133','134','135','136','137','138','139','210','211','212','213','214','215','216','217','218','219','220','221','222','223','224','225','226','227','228','229','230','231','232','233','234','235','236','237','238','239','301','302','303','304','305','306','307','401','402','403','404','405','406','502','503','504','505','507','601','602','603','604','605','606','701','702','800','801','802','803','900','901','902','903','904','910','911','912','913','914','915','916','917','918','919')
sample_IDs <- c('090','1010','1011','1012','1013','1014','1015','1016','1017','1018','1019','100','101','102','103','104','105','106','107','108','109','110','111','112','113','114','115','116','117','118','119','120','121','122','123','124','125','126','127','128','129','130','131','132','133','134','135','136','137','138','139','210','211','212','213','214','215','216','217','218','219','220','221','222','223','224','225','226','227','228','229','230','231','232','233','234','235','236','237','238','239','301','302','303','304','305','306','307','401','402','403','404','405','406','502','503','504','505','507','601','602','603','604','605','606','701','702','800','801','802','803','900','901','902','903','904','910','911','912','913','914','915','916','917','918','919')


chr01_vcf_in <- read.table(gzfile('09.nofilter.recode.nocolon.vcf.gz'), header = FALSE)
names(chr01_vcf_in) <- vcf_header

#add ancestral state to new column 
Gkirkii_seq <- read.table('roast.09.msa.Gkirkii.seqlines', header = FALSE)
chr01_vcf_in$AnState <- Gkirkii_seq[chr01_vcf_in$POS,1]
chr01_vcf_filter1 <- chr01_vcf_in[chr01_vcf_in$AnState != "-",]
chr01_vcf_filter1$AnState <- droplevels(chr01_vcf_filter1$AnState)


#add GERP score to new column (same function as before, different input?)
gerp_score <- read.table('roast.09.msa.in.gaps.rates', header = FALSE)
names(gerp_score) <- c('Neutral Rate','RS Score')
chr01_vcf_filter1$GERPgap0 <- gerp_score[chr01_vcf_filter1$POS,2]

gerp_score <- read.table('roast.09.msa.in.new.rates', header = FALSE)
names(gerp_score) <- c('Neutral Rate','RS Score')
chr01_vcf_filter1$GERP <- gerp_score[chr01_vcf_filter1$POS,2]

#determine if reference == ancestral OR alternative == ancestral
#exclude all others 
RefEqual <- chr01_vcf_filter1[c(chr01_vcf_filter1$AnState == chr01_vcf_filter1$REF),]
AltEqual <- chr01_vcf_filter1[c(chr01_vcf_filter1$AnState == chr01_vcf_filter1$ALT),]

#sub 0/0 for 2 if alternative == ancestral
AltEqual[10:135] <- lapply(AltEqual[10:135], factor, levels=c('./.', '0/0', '0/1', '1/1'), labels=c(NA, 2, 1, 0))


#sub 1/1 for  2 if reference == ancestral
RefEqual[10:135] <- lapply(RefEqual[10:135], factor, levels=c('./.', '0/0', '0/1', '1/1'), labels=c(NA, 0, 1, 2))

#recombine dfs
PostCounts <- rbind(AltEqual, RefEqual)
PostCounts <- subset(PostCounts, select=-c(QUAL, INFO, FILTER, FORMAT, ID))

#write output to file 
write.table(PostCounts, "09.DerivedCount.GERP.txt", sep='\t', row.names=F, col.names = F, quote=F)
