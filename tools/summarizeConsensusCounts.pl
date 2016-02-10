#!/usr/bin/perl
use strict;

=pod 

Summarize the output of something like:

  BOARAT/generateConsensusFromBAM.pl --bam=HifiHifi.5plxHIFIhifi_S4_L001_R1_001.fastq.gz.sorted.bam 
  
as part of a downsampling experiment.  i.e.

  BOARAT/generateConsensusFromBAM.pl --bam=HifiHifi.5plxHIFIhifi_S4_L001_R1_001.fastq.gz.sorted.bam |
		perl ../BOARAT/tools/summarizeConsensusCounts.pl --downsample=1.0

		chr3:178935995  3       6       1.0
		chr17:7578174   3       3675    1.0
		etc.



=cut

use Getopt::Long;

my %args = ();
GetOptions(
	"downsample=s" => \$args{"downsample"},
	"test" => \$args{"test"}   # test mode
);


my @results = ();
while (<>) {
	chomp;
	my ($chrom, $chromStart, $barcode, $depth, $cons, $consLong) = split /\t/;
	my $amplicon = "$chrom:$chromStart";
	
	push(@results, [ ($amplicon, $depth) ]);
}

# output some data
my @depths = qw/3 5 7 10 20 30/; 
foreach my $mindepth ( @depths ) {

	my %ampCounts = ();
	
	foreach my $r ( @results ) {
		$ampCounts{$r->[0]}++ if ( $r->[1] >= $mindepth );
		}

	foreach my $k ( keys %ampCounts ) {
		printf("%s\t%s\t%s\t%s\n", $k, $mindepth, $ampCounts{$k}, $args{"downsample"});
		}
}

exit;

