#!/usr/bin/perl
use strict;

=pod

Summarise what is in the off-target amplicons.

	time perl $BOARAT/classifyOffTargetReads.pl --fastqgz=Sample_1_ad4.R1.fastq.gz \
		--ampTable=$BOARAT/amplicon_tables/15plex_amplicons.txt --test --dump 2>&1 |
		perl reportOffTargetAmpliconComposition.pl > Sample01.OffTargetAmplicons_report.txt
		

For reference, the ideas came from my July 17, 2014 email.

Author: Paul Krzyzanowski
Email: paul.krzyzanowski@oicr.on.ca

=cut


my %f = ();
my $totalDepth = 0;
my %barcodeLocations = ();
my %primerSets = ();

while (<>) {

	chomp;

	if (/^BarcodeSuffixFound/) {
	
		$totalDepth++;
		my ($desc, $barcode, $barcodeLoc) = split /\t/;
		$barcodeLocations{$barcodeLoc}++;
	
	} elsif (/^OffTargetAmplicon/) {

		my ($desc, $amplicon, $barcodeLoc, $barcode, $primers, $seq) = split /\t/;
		$primerSets{$primers}++;
		
	}
	
		
}

print "TotalDepth: $totalDepth\n";

foreach my $bcl ( sort keys %barcodeLocations ) {
	print "BarcodeLocationCount:\t$bcl\t" . $barcodeLocations{$bcl} . "\n";
	}
	
foreach my $psl ( sort { $primerSets{$a} <=> $primerSets{$b} } keys %primerSets ) {
	print "PrimerSetCount:\t$psl\t" . $primerSets{$psl} . "\n" if ($primerSets{$psl} > 2);
	}
	


exit;

