#!/usr/bin/perl
use strict;

=pod 

Summarize per amplicon consensus depths from consensus sequences gzip file

  gunzip -c $SAMPLENAME.UIDdepths.txt.gz | perl $BHOME/tools/summarizeAmpliconConsensusDepths.pl --sampleID=$SAMPLENAME --depths=1,3,10,20,30,100 > $SAMPLENAME.consensusStatistics.txt

=cut

use Getopt::Long;
use Data::Dumper;

my %args = ();
GetOptions(
	"sampleID=s" => \$args{"sampleID"},   # sampleID
	"depths=s" => \$args{"depths"},   # sampleID
	"test" => \$args{"test"}   # test mode
);

my @depths = split(",", $args{"depths"});
my $raw_reads = '';
my %a_data = ();  # Amplicon data structure

while (<>) {

	my ($amplicon, $barcode, $barcode_count, $seq, $seqMarkup) = split /\t/;

	foreach my $d ( @depths ) {
		$a_data{$amplicon}->{$d}++ if ( $barcode_count >= $d );
	}

}

print $args{"sampleID"} . " Consensus Reads at Depth:\t" . join("\t", @depths) . "\n";
foreach my $amplicon ( keys %a_data ) {
	printf ("%s %s", $args{"sampleID"}, $amplicon);
		foreach my $d ( @depths ) {
		printf ("\t%d", $a_data{$amplicon}->{$d});
	}
	print "\n";
}

print "\n";

exit;
