#!/usr/bin/perl
use strict;

=pod

Combine HaloplexHS barcode fastqs with read fastqs

  perl combineHaloplexFastqs.pl --R1file [] --R2file [] --barcodefastq [] --outputdir []

=cut

use Getopt::Long;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use Debarcer;
use Config::General qw(ParseConfig);

my %args = ();
GetOptions(
	"outputdir=s" => \$args{"outputdir"},
	"R1file=s" => \$args{"R1file"},
	"R2file=s" => \$args{"R2file"},   
	"barcodefastq=s" => \$args{"barcodeFq"}   
);

my %haloplexBarcodes = ();
open BARCODES, "gunzip -c ".$args{"barcodeFq"}." |";
while (<BARCODES>) {
	if (/^@/) {  # If we've hit a read entry the next line is the HaloplexHS barcode
		my @readHeader = split(/\s/, $_);
		my $readID = shift @readHeader;
		my $barcode = <BARCODES>;
		chomp $barcode;
		$haloplexBarcodes{$readID} = $barcode;
	}
}
close BARCODES;
print STDERR keys(%haloplexBarcodes) . " HaloplexHS barcodes found\n";


