#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 26   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gatk"   # job name
#SBATCH --mem=900G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=GATK_filter.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel 
module load gatk/3.8-0-6lrc5l4
#module load picard/2.17.0-ft5qztz

#picard CreateSequenceDictionary R=Ghirsutum_458_v1.0.26chr.fa O=Ghirsutum_458_v1.0.26chr.dict

parallel gatk -T VariantFiltration -R Ghirsutum_458_v1.0.26chr.fa -V {}.vcf --genotypeFilterExpression "QD\<2.0\|\|FS\>60.0\|\|MQ\<40.0\|\|SOR\>4.0\|\|MQRankSum\<-12.5\|\|ReadPosRankSum\<-8.0" --genotypeFilterName "all_snp_filter" -o {}.filtered_snp.vcf "2> {}.VariantFiltration_snp.err" ::: {01..26}

parallel gatk -T SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.vcf --excludeFiltered -o {}.PASS_snp.vcf "2> {}.PASS_snp.err" ::: {01..26}

parallel gatk -T SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.vcf ----setFilteredGtToNocall -o {}.Nocall.vcf "2> {}.Nocall.err" ::: {01..26}
