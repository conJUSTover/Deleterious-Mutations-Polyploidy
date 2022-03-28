#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 30   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "bcftools"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-bcftools-merge.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load vcftools 
module load bcftools/1.9-womp5fh
module load parallel


#for i in Ghirsutum*Nocall.vcf; do vcftools --vcf $i --remove-indels --min-alleles 2 --max-alleles 2 --recode --out ${i%.Nocall.vcf}.biallele; done 
#for i in Ghirsutum*Nocall.20.0.vcf; do vcftools --vcf $i --remove-indels --min-alleles 2 --max-alleles 2 --recode --out ${i%Nocall.20.0.vcf}20.0.biallele; done 

#for i in Ghirsutum_458_v1.0.26chr.[21].[10].20.0.biallele.recode.vcf; do bcftools view -I $i -O z -o ${i%vcf}trial.vcf.gz && bcftools index ${i%vcf}trial.vcf.gz; done 

#bcftools merge -f .,PASS,ExHet20.0 -o Ghirsutum_458_v1.0.26chr.1-0.2-1.20.0.biallele.recode.vcf.gz -O z --threads 30 Ghirsutum_458_v1.0.26chr.[21].[10].20.0.biallele.recode.trial.vcf.gz

#Split Pass Names

#for i in *Nocall.[21].[01].20.0.vcf; do vcftools --vcf $i --remove-indels --min-alleles 2 --max-alleles 2 --recode --out ${i%vcf}biallele; done 
#for i in *Nocall.[21].[01].20.0.biallele.recode.vcf; do bcftools view -I $i -O z -o $i.gz && bcftools index $i.gz; done 
#bcftools merge -o separateCalls.vcf.gz -O z --threads 30 *Nocall.[21].[01].20.0.biallele.recode.vcf.gz

#for i in *NoCall.biallele.recode.vcf; do bcftools view -I $i -O z -o $i.gz && bcftools index $i.gz; done 
#bcftools merge -o separateCalls.vcf.gz -O z --threads 30 *NoCall.biallele.recode.vcf.gz

#for i in {01..26}; do vcftools --gzvcf separateCalls.vcf.gz --out $i.nofilter --chr $i --recode; done 
#parallel 'vcftools --gzvcf separateCalls.vcf.gz --out {}.nofilter --chr {} --recode' ::: {01..26}
 
ls *0.NoCall.biallele.recode.vcf.gz | sed 's/.NoCall.biallele.recode.vcf.gz//g'| parallel 'vcftools --gzvcf {}.NoCall.biallele.recode.vcf.gz --out {} --freq2'

