#!/usr/bin/bash 

module load vcftools/0.1.14-xp36ajy
module load bcftools/1.9-womp5fh
module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
module load parallel

parallel gatk VariantFiltration -R Ghirsutum_458_v1.0.26chr.fa -V {}.prelim.recode.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --output {}.filtered_all.vcf "2> {}.VariantFiltration_all.err" ::: {01..26}

parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_all.vcf --set-filtered-gt-to-nocall --output {}.all.GATKfiltered.vcf 2> {}.all.GATKfiltered.er ::: {01..26}

bcftools concat -o all.GATK.filtered.vcf *all.GATKfiltered.vcf

