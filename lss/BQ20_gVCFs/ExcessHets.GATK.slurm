#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 36:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 30   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "AD2 gatk4"   # job name
#SBATCH --mem=90G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-GATK4-21.ExHet.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,singularity,gpu #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel
module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
module load vcftools/0.1.14-xp36ajy
module load picard/2.17.0-ft5qztz
module load sentieon-genomics/201808.01-opfuvzr

ref="Ghirsutum_458_v1.0.26chr"
pop="2.1"

#sentieon driver -t 30 -r $ref.fa --algo GVCFtyper $ref.$pop.vcf 21[0-9].gVCF &> $pop.vcf.timelog

#ls $ref.fa | parallel gatk VariantFiltration -R {} -V {.}.$pop.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --output {.}.$pop.filtered_snp.vcf "2> {.}.$pop.VariantFiltration_snp.err"

#ls $ref.fa | parallel gatk SelectVariants -R {} --variant {.}.$pop.filtered_snp.vcf --set-filtered-gt-to-nocall --output {.}.$pop.Nocall.vcf "2> {.}.$pop.Nocall.err"

#ls $ref.fa | parallel vcftools --vcf {.}.$pop.Nocall.vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --out {.}.$pop.biallele

#echo "First Pass Done!" 


ExHet=20.0

#ls $ref.fa | parallel gatk VariantFiltration -R {} -V {.}.$pop.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --filter-expression "ExcessHet\>$ExHet" --filter-name "ExHet$ExHet" --output {.}.$pop.filtered_snp.$ExHet.vcf "2> {.}.$pop.VariantFiltration_snp.$ExHet.err"

#ls $ref.fa | parallel gatk SelectVariants -R {} --variant {.}.$pop.filtered_snp.$ExHet.vcf --set-filtered-gt-to-nocall --output {.}.$pop.Nocall.$ExHet.vcf "2> {.}.$pop.Nocall.$ExHet.err"

#ls $ref.fa | parallel vcftools --vcf {.}.$pop.Nocall.$ExHet.vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --out {.}.$pop.$ExHet.biallele

#for i in $ref.*biallele.recode.vcf; do vcftools --vcf $i --remove-filtered-all --min-alleles 2 --max-alleles 2 --out ${i%.vcf} --freq2; vcftools --vcf $i --remove-filtered-all --min-alleles 2 --max-alleles 2 --out ${i%.vcf} --window-pi 50000; vcftools --vcf $i --remove-filtered-all --min-alleles 2 --max-alleles 2 --out ${i%.vcf} --site-pi;done 




ls $ref.fa | parallel gatk VariantFiltration -R {} -V {.}.$pop.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --filter-expression "ExcessHet\>$ExHet" --filter-name "ExHet$pop.$ExHet" --output {.}.$pop.filtered_snp.$pop.$ExHet.vcf "2> {.}.$pop.VariantFiltration_snp.$pop.$ExHet.err"

ls $ref.fa | parallel gatk SelectVariants -R {} --variant {.}.$pop.filtered_snp.$pop.$ExHet.vcf --set-filtered-gt-to-nocall --output {.}.$pop.Nocall.$pop.$ExHet.vcf "2> {.}.$pop.Nocall.$pop.$ExHet.err"

pop="1.0"

ls $ref.fa | parallel gatk VariantFiltration -R {} -V {.}.$pop.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --filter-expression "ExcessHet\>$ExHet" --filter-name "ExHet$pop.$ExHet" --output {.}.$pop.filtered_snp.$pop.$ExHet.vcf "2> {.}.$pop.VariantFiltration_snp.$pop.$ExHet.err"

ls $ref.fa | parallel gatk SelectVariants -R {} --variant {.}.$pop.filtered_snp.$pop.$ExHet.vcf --set-filtered-gt-to-nocall --output {.}.$pop.Nocall.$pop.$ExHet.vcf "2> {.}.$pop.Nocall.$pop.$ExHet.err"

echo "refiltering complete"


