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
	"byAmplicon" => \$args{"byAmplicon"},
	
);

if ( $args{"byAmplicon"} ) {

	open RCFILE, $args{"consSeqFile"} || die "Need to provide a consensusSequences.txt file\n\n";
	my @consenses = <RCFILE>;
	close RCFILE;

	my %nuc_idx = ();

	foreach my $line ( @consenses ) { 
		chomp;
		my ( $amp, $bc, $depth, $readC, $readMarkup) = split /\t/, $line;
		# print "$amp $bc $depth\n";
		
		for ( my $i = 0; $i < length($bc); $i++ ) {
			my $nuc = substr($bc, $i, 1);
			
			if ($depth < 10) {
				$nuc_idx{$amp}{$i}->{$nuc}->{"low"}++;
				$nuc_idx{$amp}{$i}->{"Total"}->{"low"}++;
				} elsif ( $depth >= 100 ) {
				$nuc_idx{$amp}{$i}->{$nuc}->{"high"}++;
				$nuc_idx{$amp}{$i}->{"Total"}->{"high"}++;
				}
			
			
		}

	}

	foreach my $amp ( sort {$a cmp $b} keys %nuc_idx ) {
		
		print "# $amp Low abundance barcodes\n";
		print join("\t", qw/Amplicon Abundance Site A nA P(A) C nC P(C) G nG P(G) T nT P(T)/) . "\n";
		foreach my $i ( sort {$a <=> $b} keys %{$nuc_idx{$amp}} ) {
			printf( "%s\t%s\t%s", $amp, "Low", $i);
			foreach my $nuc ( sort {$a cmp $b} keys %{$nuc_idx{$amp}{$i}} ) {
				printf( "\t%s\t%s\t%0.3f", $nuc, $nuc_idx{$amp}{$i}->{$nuc}{"low"}, $nuc_idx{$amp}{$i}->{$nuc}{"low"} / $nuc_idx{$amp}{$i}->{"Total"}{"low"});
			}
			print "\n";
		}

		print "\n";
			
		print "# $amp High abundance barcodes\n";
		print join("\t", qw/Amplicon Abundance Site A nA P(A) C nC P(C) G nG P(G) T nT P(T)/) . "\n";
		foreach my $i ( sort {$a <=> $b} keys %{$nuc_idx{$amp}} ) {
			printf( "%s\t%s\t%s", $amp, "High", $i);
			foreach my $nuc ( sort {$a cmp $b} keys %{$nuc_idx{$amp}{$i}} ) {
				printf( "\t%s\t%s\t%0.3f", $nuc, $nuc_idx{$amp}{$i}->{$nuc}{"high"}, $nuc_idx{$amp}{$i}->{$nuc}{"high"} / $nuc_idx{$amp}{$i}->{"Total"}{"high"});
			}
			print "\n";
		}
		
		print "\n\n";
	}

} else {

open RCFILE, $args{"consSeqFile"} || die "Need to provide a consensusSequences.txt file\n\n";
my @consenses = <RCFILE>;
close RCFILE;

my %nuc_idx = ();

foreach my $line ( @consenses ) { 
	chomp;
	my ( $amp, $bc, $depth, $readC, $readMarkup) = split /\t/, $line;
	# print "$amp $bc $depth\n";
	
	for ( my $i = 0; $i < length($bc); $i++ ) {
		my $nuc = substr($bc, $i, 1);
		
		if ($depth < 10) {
			$nuc_idx{$i}->{$nuc}->{"low"}++;
			$nuc_idx{$i}->{"Total"}->{"low"}++;
			} elsif ( $depth >= 100 ) {
			$nuc_idx{$i}->{$nuc}->{"high"}++;
			$nuc_idx{$i}->{"Total"}->{"high"}++;
			}
		
		
	}

}

print "Low abundance barcodes\n";
print join("\t", qw/Site A nA P(A) C nC P(C) G nG P(G) T nT P(T)/) . "\n";
foreach my $i ( sort {$a <=> $b} keys %nuc_idx ) {
	printf( "%s", $i);
	foreach my $nuc ( sort {$a cmp $b} keys %{$nuc_idx{$i}} ) {
		printf( "\t%s\t%s\t%0.3f", $nuc, $nuc_idx{$i}->{$nuc}{"low"}, $nuc_idx{$i}->{$nuc}{"low"} / $nuc_idx{$i}->{"Total"}{"low"});
	}
	print "\n";
}

print "\n";
	
print "High abundance barcodes\n";
print join("\t", qw/Site A nA P(A) C nC P(C) G nG P(G) T nT P(T)/) . "\n";
foreach my $i ( sort {$a <=> $b} keys %nuc_idx ) {
	printf( "%s", $i);
	foreach my $nuc ( sort {$a cmp $b} keys %{$nuc_idx{$i}} ) {
		printf( "\t%s\t%s\t%0.3f", $nuc, $nuc_idx{$i}->{$nuc}{"high"}, $nuc_idx{$i}->{$nuc}{"high"} / $nuc_idx{$i}->{"Total"}{"high"});
	}
	print "\n";
}

}	
	


exit;

