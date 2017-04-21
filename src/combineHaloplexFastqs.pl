#!/usr/bin/perl
use strict;

=pod

Combine HaloplexHS barcode fastqs with read fastqs

  perl combineHaloplexFastqs.pl --R1file [] --R2file [] --barcodefastq [] --outputdir []

=cut

use Getopt::Long;
use Data::Dumper;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/";
use Debarcer;
use Config::General qw(ParseConfig);

print STDERR "Starting combineHaloplexFastqs.pl\n";

my %args = ();
GetOptions(
	"outputdir=s" => \$args{"outputdir"},
	"R1file=s" => \$args{"R1file"},
	"R2file=s" => \$args{"R2file"},   
	"barcodefastq=s" => \$args{"barcodeFq"}   
);
my %haloplexBarcodes = ();
my $outputfile = '';

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


# Add the barcodes to the Read 1 file
open INFILE, "gunzip -c ".$args{"R1file"}." |";
$outputfile = $args{"outputdir"}."/". basename($args{"R1file"});
$outputfile =~ s/_R\d_/_R1_/;  # Ensure the file indicated R1
$outputfile =~ s/\.gz$//; # output uncompressed fastq
open OUTFILE, ">$outputfile";
while (<INFILE>) {
	if (/^@/) {  # If we've hit a read entry the next line is the HaloplexHS barcode
		my @readHeader = split(/\s/, $_);
		$readHeader[0] .= ":HaloplexHS-".$haloplexBarcodes{$readHeader[0]};
		print OUTFILE join(" ", @readHeader) . "\n";
	} else {
		print OUTFILE;
	}
}
close INFILE;
close OUTFILE;
system("gzip $outputfile");
print STDERR "Added HaloplexHS barcodes to R1 file\n";
 