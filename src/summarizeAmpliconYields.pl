#!/usr/bin/perl
use strict;

=pod

Summarize amplicon yields from information in the DeBarcEr.log file.

  cat DeBarcEr.log | perl summarizeAmpliconYields.pl

=cut

use Getopt::Long;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use Debarcer;
use Config::General qw(ParseConfig);

my %args = ();
GetOptions(
	"config=s" => \$args{"configfile"},
	"sampleID=s" => \$args{"sampleID"},   # sampleID
	"test" => \$args{"test"}   # test mode
);

# Section to load parameters from a Config::Simple format file [Not used here yet. FIXME]
die "Need to supply a config file.\n" unless ( $args{"configfile"} );
my %config = ();
%config = ParseConfig($args{"configfile"});
my $ampliconTable = ( $config{"ampliconTable"} ) ? $config{"ampliconTable"} : "$FindBin::Bin/../data/all_amplicons.txt" ;

my $raw_reads = '';
my %a_data = ();  # Amplicon data structure

my %ampliconInfo = ();
&Debarcer::loadAmpliconData($ampliconTable, \%ampliconInfo); 
my %positionAliases = &Debarcer::getPositionAliases($ampliconTable);

while (<>) {

	if (/Raw reads mapped by bwa/) {
		$_ =~ m/(\d+)$/g;
		$raw_reads = $1;
		}

	if (/depth\|count\|coverage/) {
		chomp;
		my ($amplicon, $text, $depth, $count, $coverage) = split /\t/;
		$a_data{$amplicon}->{$depth}{"count"} = $count;
		$a_data{$amplicon}->{$depth}{"coverage"} = $coverage;
		$a_data{"Aggregate"}->{$depth}{"count"} += $count;
		$a_data{"Aggregate"}->{$depth}{"coverage"} += $coverage;
	
	}

}

print "Debarcer Summary Statistics for:\t" . $args{"sampleID"} . "\n\n";

print "Raw reads:\t$raw_reads\n";
printf ("Amplicon Reads Identified:\t%d\n", $a_data{"Aggregate"}->{1}{"coverage"});
printf ("Amplicon Yield:\t%0.6f\n", $a_data{"Aggregate"}->{1}{"coverage"} / $raw_reads);
print "\n";

foreach my $d ( qw/1 3 10 20 30 100/ ) {
	printf ("Consensus(%d) Reads:\t%s\n", $d, ( exists $a_data{"Aggregate"}->{$d}{"count"} ) ? $a_data{"Aggregate"}->{$d}{"count"} : "NA" );
 }
print "\n";

foreach my $d ( qw/1 3 10 20 30 100/ ) {
	printf ("Overall Consensus(%d) Yield:\t%s\n", $d, 
		( exists $a_data{"Aggregate"}->{$d}{"count"} ) ? sprintf("%0.6f", $a_data{"Aggregate"}->{$d}{"count"} / $raw_reads) : "NA");
	printf ("Consensus(%d) Yield:\t%s\n", $d, 
		( exists $a_data{"Aggregate"}->{$d}{"count"} ) ? sprintf("%0.6f", $a_data{"Aggregate"}->{$d}{"count"} / $a_data{"Aggregate"}->{1}{"coverage"} ) : "NA");
}
print "\n";

print "All amplicons identified:\n";
foreach my $a ( sort keys %a_data ) {
	my $ampliconName = ( exists $positionAliases{$a} ) ? $positionAliases{$a} : "-";
	printf( "%s\t%s\t%d\n", $a, $ampliconName, $a_data{$a}->{1}{"coverage"});
}

exit;

