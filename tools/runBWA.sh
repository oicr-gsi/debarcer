#!/bin/bash

. /etc/profile
module load samtools/0.1.19;
module load bwa/0.7.12;

FASTQ=$1
SAMPLENAME=$2
echo "Running bwa/0.7.12 on Sample $SAMPLENAME using fastq: $FASTQ";

HG19="/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.6.2/hg19_random.fa"
NPROCS=`grep processor /proc/cpuinfo | wc -l`
echo "bwa using $NPROCS processors";

bwa mem -t $NPROCS $HG19 $FASTQ | samtools view -bS - | samtools sort - $SAMPLENAME.$FASTQ.sorted
samtools index $SAMPLENAME.$FASTQ.sorted.bam
# bwa mem $HG19 $FASTQ > $SAMPLENAME.$FASTQ.sam	

