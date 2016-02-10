#!/usr/bin/perl
use strict;

=pod

Summarise the numbers of analyzable amplicon-barcode groups
at a range of cutoffs.

	time gunzip -c $SAMPLENAME.UIDdepths.txt.gz |
		perl reportBarcodeComposition.pl > Sample01.barcode_report.txt

For reference, the ideas came from my July 17, 2014 email.

Author: Paul Krzyzanowski
Email: paul.krzyzanowski@oicr.on.ca

=cut



my  %f = ();
my $totalDepth = 0;

while (<>) {

	chomp;
	my ($amplicon, $barcode, $depth) = split(/\t/);
	my $family = $amplicon . "_" . $barcode;
	$f{$family} = $depth;
	$totalDepth += $depth;
	
}

my @families = keys %f;
my $totalFamilies = scalar @families;

print join("\t", qw/MinDepth FamiliesAtDepth FamiliesAtDepthOrGreater ReadsAtDepth PercReadsAtDepth CumulativeReads CumulativePercOfTotalReads/) . "\n";
foreach my $d ( 1..500 ) {

	my $FamiliesWithDepth = 0;
	my $TotalFamiliesWithMinDepth = 0;
	my $CumulativeReads = 0;
	foreach my $fam ( keys %f ) { 
		$FamiliesWithDepth++ if ( $f{$fam} == $d ); 
		$TotalFamiliesWithMinDepth++ if ( $f{$fam} >= $d ); 
		$CumulativeReads += $f{$fam} if ( $f{$fam} >= $d ); 
		}
	my $Reads_at_d = $d * $FamiliesWithDepth;
	my $PercReadsAtDepth = $Reads_at_d / $totalDepth;
	my $CumReadsAsPercOfTotal = $CumulativeReads / $totalDepth;
	
	printf("%s\t%s\t%s\t%s\t%0.6f\t%s\t%0.6f\n", $d, $FamiliesWithDepth, $TotalFamiliesWithMinDepth, $Reads_at_d, $PercReadsAtDepth, $CumulativeReads, $CumReadsAsPercOfTotal);
	
}



exit;


# SVNID: $Id: reportBarcodeComposition.pl 324 2015-06-10 18:10:33Z pkrzyzanowski $
