#!/usr/bin/perl
use strict;

=pod

This script reads in family barcodes and their depths and derives information 
on the kmer composition identified therein.

(c) Paul Krzyzanowski 2015
Email: paul.krzyzanowski@oicr.on.ca

=cut

use List::Util qw/max shuffle/;
use Getopt::Long;

my %args = ();
GetOptions(
	"infile=s" => \$args{"consSeqFile"},
	"byAmplicon" => \$args{"byAmplicon"}
	
);

if ( $args{"byAmplicon"} ) {

	open RCFILE, $args{"consSeqFile"} || die "Need to provide a consensusSequences.txt file\n\n";
	my @consenses = <RCFILE>;
	close RCFILE;

	my %kmer_idx = ();

	foreach my $line ( @consenses ) { 
		chomp;
		my ( $amp, $bc, $depth, $readC, $readMarkup) = split /\t/, $line;
		# print "$amp $bc $depth\n";
		
		for ( my $i = 0; $i <= (length($bc) - 3); $i++ ) {
			my $kmer = substr($bc, $i, 3);
			
			# $kmer_idx{$kmer}->{$depth}++;
			
			if ($depth < 10) {
				$kmer_idx{$amp}{$kmer}->{"low"}++;
				$kmer_idx{$amp}{"Total"}->{"low"}++;
				} elsif ( $depth >= 100 ) {
				$kmer_idx{$amp}{$kmer}->{"high"}++;
				$kmer_idx{$amp}{"Total"}->{"high"}++;
				}
			
			
		}

	}

	foreach my $amp ( sort {$a cmp $b}  keys %kmer_idx ) {
	foreach my $kmer ( sort {$a cmp $b} keys %{$kmer_idx{$amp}} ) {
		my $nLow = $kmer_idx{$amp}{$kmer}->{"low"};
		my $nHigh = $kmer_idx{$amp}{$kmer}->{"high"};
		my $pLow = $nLow / $kmer_idx{$amp}{"Total"}->{"low"};
		my $pHigh = $nHigh / $kmer_idx{$amp}{"Total"}->{"high"};
		printf( "%s\t%s\t%s\t%0.4f\t%s\t%0.4f\t%0.4f\n", $amp, $kmer, $nLow, $pLow, $nHigh, $pHigh, $pHigh - $pLow);
		}
	}

} else {
	
open RCFILE, $args{"consSeqFile"} || die "Need to provide a consensusSequences.txt file\n\n";
my @consenses = <RCFILE>;
close RCFILE;

my %kmer_idx = ();

foreach my $line ( @consenses ) { 
	chomp;
	my ( $amp, $bc, $depth, $readC, $readMarkup) = split /\t/, $line;
	# print "$amp $bc $depth\n";
	
	for ( my $i = 0; $i <= (length($bc) - 3); $i++ ) {
		my $kmer = substr($bc, $i, 3);
		
		# $kmer_idx{$kmer}->{$depth}++;
		
		if ($depth < 10) {
			$kmer_idx{$kmer}->{"low"}++;
			$kmer_idx{"Total"}->{"low"}++;
			} elsif ( $depth >= 100 ) {
			$kmer_idx{$kmer}->{"high"}++;
			$kmer_idx{"Total"}->{"high"}++;
			}
		
		
	}

}

foreach my $kmer ( sort {$a cmp $b} keys %kmer_idx ) {
	my $nLow = $kmer_idx{$kmer}->{"low"};
	my $nHigh = $kmer_idx{$kmer}->{"high"};
	my $pLow = $nLow / $kmer_idx{"Total"}->{"low"};
	my $pHigh = $nHigh / $kmer_idx{"Total"}->{"high"};
	printf( "%s\t%s\t%0.4f\t%s\t%0.4f\t%0.4f\n", $kmer, $nLow, $pLow, $nHigh, $pHigh, $pHigh - $pLow);
	}
	
}


exit;

