# use this simple file as Sentieon.sh base_name_of_sample number_of_threads reference

module load bwa/0.7.17-zhcbtza

# generate this index before loading Sentieon if you ever want the index to finish
# time ( bwa index A2Du_26.fasta )

#module use /opt/rit/local-test/modules
module load sentieon-genomics/201808.01-opfuvzr
module load samtools/1.9-k6deoga
#module load trimmomatic/0.36-lkktrba

# samtools faidx A2Du_26.fasta

# so this command doesn't work
# it doesn't like the || echo part
# it also doesn't like to pipe directly into sentieon
# (bwa mem -M -K 10000000 -R "@RG\tID:$1 \tSM:$1 \tPL:ILLUMINA" -t $2 $3 $1.1.fq.gz $1.2.fq.gz || echo -n 'error' )| sentieon util sort -o $1.sort.bam -t $2 --sam2bam -i

#time(trimmomatic PE -threads $2 $1.1.fq.gz $1.2.fq.gz $1.1.paired.fq.gz $1.1.unpaired.fq.gz $1.2.paired.fq.gz $1.2.unpaired.fq.gz ILLUMINACLIP:Adapters.fa:2:30:15 LEADING:28 TRAILING:28 SLIDINGWINDOW:8:28 SLIDINGWINDOW:1:10 MINLEN:65 TOPHRED33) &> $1.trimmomatic.timelog

time ( bwa mem -M -K 10000000 -R "@RG\tID:$1 \tSM:$1 \tPL:ILLUMINA" -t $2 $3 $1.1.paired.fq.gz $1.2.paired.fq.gz > $1.sam  ) &> $1.bwa.timelog

samtools view -F 260 -bu -@ $2 $1.sam | samtools sort -m 4G -@ $2 -o $1.sort.bam && rm $1.sam
samtools index -@ $2 $1.sort.bam

#sorting done in previous step; not needed here
#time ( sentieon util sort -o $1.sort.bam -t $2 --sam2bam -i $1.sam  ) &> $1.sam2bam.timelog

# calculates GC bias; compute if you want
time ( sentieon driver -t $2 -r $3 -i $1.sort.bam --algo GCBias --summary $1.GC.summary $1.GC.metric --algo MeanQualityByCycle $1.MQ.metric --algo QualDistribution $1.QD.metric --algo InsertSizeMetricAlgo $1.IS.metric --algo AlignmentStat A2_060.ALN.metric ) &> $1.metric.timelog

time ( sentieon plot metrics -o $1.metric.pdf gc=$1.GC.metric mq=$1.MQ.metric qd=$1.QD.metric isize=$1.IS.metric ) &> $1.metricPlot.timelog

#### these are for removing duplicate reads
time ( sentieon driver -t $2 -i $1.sort.bam --algo LocusCollector --fun score_info $1.score ) &> $1.locusCollect.timelog

time ( sentieon driver -t $2 -i $1.sort.bam --algo Dedup --rmdup --score_info $1.score --metrics $1.dedup.metric $1.dedup.bam ) &> $1.rmDup.timelog

#### indel realigner

time ( sentieon driver -t $2 -r $3 -i $1.dedup.bam --algo Realigner $1.realign.bam ) &> $1.realign.timelog

### same as HaplotypeCaller
# right now, set only to emit SNPs
time ( sentieon driver -t $2 -r $3 -i $1.realign.bam --algo Haplotyper $1.gVCF --emit_mode gvcf ) &> $1.gvcf.timelog

# time ( sentieon driver -t 96 -r A2Du_26.fasta --algo GVCFtyper A2Du.vcf *.gVCF ) &> A2Du.vcf.timelog
# module load snphylo
# sed 's/Chr//g' A2Du.vcf > A2Du.numOnly.vcf
# snphylo.sh -v A2Du.numOnly.vcf -c 5 -P A2phylo -b -B 10000




