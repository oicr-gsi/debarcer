#!/bin/bash

echo 'Sourcing ~/.debarcer'
source ~/.debarcer

SAMBIN="$SAMTOOLSROOT/bin"
BWABIN="$BWAROOT"

FASTQ=$1
SAMPLENAME=$2
echo "Running bwa/0.7.12 on Sample $SAMPLENAME using fastq: $FASTQ";

echo "bwa using $NPROCS processors";

$BWAROOT/bwa mem -t $NPROCS $HG19 $FASTQ | 
	$SAMTOOLSROOT/bin/samtools view -bS - | 
	$SAMTOOLSROOT/bin/samtools sort - $SAMPLENAME.$FASTQ.sorted
$SAMTOOLSROOT/bin/samtools index $SAMPLENAME.$FASTQ.sorted.bam
# bwa mem $HG19 $FASTQ > $SAMPLENAME.$FASTQ.sam	

