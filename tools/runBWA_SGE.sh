#!/bin/bash
#$ -N "HFHFbwa" 
#$ -l h_vmem=32G 
#$ -cwd 



. /etc/profile
# module load picard/1.90;
 module load samtools/0.1.19;

FASTQ="5plxHIFIhifi_S4_L001_R1_001.fastq.gz"
SAMPLENAME="HifiHifi"
HG19="/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.6.2/hg19_random.fa"

module load bwa/0.7.12;
bwa mem -t `grep processor /proc/cpuinfo | wc -l` $HG19 $FASTQ | samtools view -bS > $FASTQ.bam
samtools index $FASTQ.bam



