#!/usr/bin/bash 

module load vcftools/0.1.14-xp36ajy
module load bcftools/1.9-womp5fh
module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
module load parallel 

sed "s/nan/NaN/g" $1.vcf > $1.new.vcf && echo "$1 Sed Done"

bcftools reheader --samples header.txt -o $1.vcf $1.new.vcf && rm $1.new.vcf

vcftools --vcf $1.vcf --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 --recode --out $1.prelim

echo $1 | parallel gatk VariantFiltration -R Ghirsutum_458_v1.0.26chr.fa -V {}.prelim.recode.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --output {}.filtered_snp.vcf "2> {}.VariantFiltration_snp.err"

echo "GATK variants identified"

module load python/3.8.0-k7w5uj4
python ExcessHets_bypop.py $1.filtered_snp.vcf pops.lsts > $1.filtered_snp_exhet.vcf
gatk IndexFeatureFile -F $1.filtered_snp_exhet.vcf
gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant $1.filtered_snp_exhet.vcf --set-filtered-gt-to-nocall --output $1.biallele.recode.vcf 2> $1.Nocall.err

echo "GATK filtering completed"

# Find Sites that are present in 2+ outgroups and not segregating
vcftools --vcf $1.biallele.recode.vcf --keep out.lst --freq2 --out $1.outgroup
grep -v "CHROM" $1.outgroup.frq | awk 'OFS="\t" {if($4 > 3 && $6 == 1) {print $1,$2}}' > $1.keep.pos
grep -v "CHROM" $1.outgroup.frq | awk 'OFS="\t" {if($4 > 3 && $6 == 0) {print $1,$2}}' >> $1.keep.pos

# Remove positions 
vcftools --vcf $1.biallele.recode.vcf --positions $1.keep.pos --recode --recode-INFO-all --out $1.outgroup

# Find sites that are segregating within species of interest
vcftools --vcf $1.outgroup.recode.vcf --remove out.lst --freq2 --out $1.NotSeg
grep -v "CHROM" $1.NotSeg.frq | awk 'OFS="\t" {if($6 < 1 && $6 > 0) {print $1,$2}}' > $1.NotSeg.pos

# Keep those postiions 
vcftools --vcf $1.outgroup.recode.vcf --positions $1.NotSeg.pos --recode --recode-INFO-all --out $1.ancestral

# Keep only GT information 
bcftools annotate -x ^FORMAT/GT -o $1.GT.vcf $1.ancestral.recode.vcf


cat $1.GT.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > $1.GT_sorted.vcf



