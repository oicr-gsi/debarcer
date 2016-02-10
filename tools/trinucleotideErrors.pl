#!/usr/bin/perl
use strict;

=pod

This script reads rawPositionalComposition.txt files and derives information 
on the error profile identified therein.

(c) Paul Krzyzanowski 2014
Email: paul.krzyzanowski@oicr.on.ca

=cut

use List::Util qw/max shuffle/;
use Getopt::Long;

my %args = ();
GetOptions(
	"infile=s" => \$args{"rawCountFile"}
);

open RCFILE, $args{"rawCountFile"} || die "Need to provide a PositionalComposition.txt file\n\n";
my @rawCounts = <RCFILE>;
close RCFILE;

my %positionBases = ();
&populatePositions(\%positionBases, \@rawCounts);

# &outputSNPrates(\%positionBases, \@rawCounts);

&outputTrinucRates(\%positionBases, \@rawCounts);


exit;



#########################################################################################

sub populatePositions {
	my ($href, $aref) = @_;
	
	foreach ( @$aref ) {
		chomp;
		my ($amp, $loc, $a, $c, $t, $g, $depth) = split /\t/;
		my $ampPosition = "$amp-$loc";
	
		# decide which base it is from the read data
		my $maxbase = max ($a, $c, $t, $g);
		my $base = '';
		if ( $a == $maxbase ) {
			$base = "A";
		} elsif ( $c == $maxbase ) {
			$base = "C";
		} elsif ( $t == $maxbase ) {
			$base = "T";
		} elsif ( $g == $maxbase ) {
			$base = "G";
		} else {
			die "Can't decide on a base for $amp position $loc\n";
		}
		
		$href->{$amp}{$loc} = $base;
	
	}
	
}
	

sub outputSNPrates {
	my ($href, $aref) = @_;
	
	my %SNPrates = ();
	
	foreach ( @$aref ) {
		chomp;
		my ($amp, $loc, $a, $c, $t, $g, $depth) = split /\t/;
		my $refBase = $href->{$amp}{$loc};
		
		$SNPrates{$refBase}->{"A"} += $a;
		$SNPrates{$refBase}->{"C"} += $c;
		$SNPrates{$refBase}->{"T"} += $t;
		$SNPrates{$refBase}->{"G"} += $g;
		$SNPrates{$refBase}->{"Depth"} += $depth;
		$SNPrates{$refBase}->{"Count"}++;
		
	}
	
	
# output some small tables.
	# print "# Raw counts\n";
	# print join("\t", qw/refBase A C T G totalDepth nObserved/) . "\n";
	# foreach my $refBase ( qw/A C T G/ ) {
		# print "$refBase";
		# foreach my $nuc ( qw/A C T G/ ) {
			# print "\t" . $SNPrates{$refBase}->{$nuc};
		# }
		# print "\t" . $SNPrates{$refBase}->{"Depth"};
		# print "\t" . $SNPrates{$refBase}->{"Count"};
		# print "\n";
	# }

	print "# Normalized counts, per million\n";
	print join("\t", qw/refBase A C T G totalDepth nObserved/) . "\n";
	foreach my $refBase ( qw/A C T G/ ) {
		print "$refBase";
		foreach my $nuc ( qw/A C T G/ ) {
			printf("\t%0.0f", $SNPrates{$refBase}->{$nuc} / $SNPrates{$refBase}->{"Depth"} * 1000000) ;
		}
		print "\t" . $SNPrates{$refBase}->{"Depth"};
		print "\t" . $SNPrates{$refBase}->{"Count"};
		print "\n";
	}

	
}
	

sub outputTrinucRates {
	my ($href, $aref) = @_;
	
	my %SNPrates = ();
	
	foreach ( @$aref ) {
		chomp;
		my ($amp, $loc, $a, $c, $t, $g, $depth) = split /\t/;
		my $refBase = $href->{$amp}{$loc};
		
		# Skip if we can't make a trinucleotide at this position
		next unless ( exists $href->{$amp}{($loc-1)} );
		next unless ( exists $href->{$amp}{($loc+1)} );
		
		my $upBase = $href->{$amp}{($loc-1)};
		my $downBase = $href->{$amp}{($loc+1)};
		my $refTrinuc = $upBase . $refBase . $downBase;
				
		$SNPrates{$refTrinuc}->{"A"} += $a;
		$SNPrates{$refTrinuc}->{"C"} += $c;
		$SNPrates{$refTrinuc}->{"T"} += $t;
		$SNPrates{$refTrinuc}->{"G"} += $g;
		$SNPrates{$refTrinuc}->{"Depth"} += $depth;
		$SNPrates{$refTrinuc}->{"Count"}++;
		
	}
	
	
# output some small tables.
	# print "# Raw counts\n";
	# print join("\t", qw/refBase A C T G totalDepth nObserved/) . "\n";
	# foreach my $refBase ( sort keys %SNPrates ) {
		# print "$refBase";
		# foreach my $nuc ( qw/A C T G/ ) {
			# print "\t" . $SNPrates{$refBase}->{$nuc};
		# }
		# print "\t" . $SNPrates{$refBase}->{"Depth"};
		# print "\t" . $SNPrates{$refBase}->{"Count"};
		# print "\n";
	# }

	print "# Normalized counts, per million\n";
	print join("\t", qw/refBase A C T G totalDepth nObserved/) . "\n";
	foreach my $refBase ( sort keys %SNPrates ) {
		print "$refBase";
		foreach my $nuc ( qw/A C T G/ ) {
			printf("\t%0.0f", $SNPrates{$refBase}->{$nuc} / $SNPrates{$refBase}->{"Depth"} * 1000000) ;
		}
		print "\t" . $SNPrates{$refBase}->{"Depth"};
		print "\t" . $SNPrates{$refBase}->{"Count"};
		print "\n";
	}

	
}
	

# SVNID: $Id$
