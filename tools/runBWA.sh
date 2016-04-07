#!/bin/bash

echo 'Sourcing ~/.debarcer'
source ~/.debarcer

FASTQGZ=$1
SAMPLENAME=$2
FASTQBASENAME=${FASTQGZ##*/}
SAMPLEPREFIX=$SAMPLENAME"."$FASTQBASENAME

echo "Running bwa/0.7.12 on Sample $SAMPLENAME using fastq: $FASTQGZ";
echo "bwa using $NPROCS processors";

$BWAROOT/bwa mem -t $NPROCS $HG19 $FASTQGZ | 
	$SAMTOOLSROOT/bin/samtools view -bS - | 
	$SAMTOOLSROOT/bin/samtools sort - $SAMPLEPREFIX.sorted
$SAMTOOLSROOT/bin/samtools index $SAMPLEPREFIX.sorted.bam
# bwa mem $HG19 $FASTQ > $SAMPLENAME.$FASTQ.sam	

