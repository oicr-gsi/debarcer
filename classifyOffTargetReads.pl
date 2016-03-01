#!/usr/bin/perl
use strict;

=pod

=head1 classifyOffTargetReads.pl

Description

=head2 Usage

   perl classifyOffTargetReads.pl --dump --fastqgz=../fastqs/Sample_01.R1.fastq.gz --ampTable=$AMPLICON_TABLE
   	   
=head2 Author

Paul M Krzyzanowski
pmkrzyzanowski@gmail.com
(c) 2014-2015

=cut


# Setup general modules
# use Cwd 'abs_path';
# use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use Debarcer;
use Getopt::Long;

my %args = ();

GetOptions(
	"fastqgz=s" => \$args{"fastqgz"},
	"ampTable=s" => \$args{"ampTable"},  # Amplicon table filename
	"dump" => \$args{"dump"},  # dump mode
	"test" => \$args{"test"}   # test mode
);

# Put the amplicon table into a global array
unless ( $args{"ampTable"} ) { # If no argument given use this default file
	$args{"ampTable"} =  '/u/pkrzyzanowski/projects/EAC/molecular_barcoding/svn/trunk/debarcer/amplicon_tables/2014July_5plex_amplicons.txt'; 
	}
open AMPS, $args{"ampTable"} or die "Can't open table of amplicon information\n";
my @ampliconsFile = <AMPS>;
# print STDERR "Amplicons:\n@ampliconsFile\n";
close AMPS;

%Debarcer::primerSets = &Debarcer::listPrimers(\@ampliconsFile);
%Debarcer::ampSeq = &Debarcer::ampliconSequences(\@ampliconsFile);
%Debarcer::ampLen = &Debarcer::ampliconLengths(\@ampliconsFile);

my %readData = ();


my $infile = $args{"fastqgz"};
if ( $args{"test"} ) {
	open INFILE, "gunzip -c $infile | head -n 400000 |";
	} else {
	open INFILE, "gunzip -c $infile |";
	}

if ( $args{"dump"} ) { # Run dump mode

	while (my $header = <INFILE>) {
		my $seq = <INFILE> || die; 
		next if ($header =~ /\+/); # Skip since $seq is the string of Q-values.
		my $separator = <INFILE> || die; # Read the + symbol
		my $quals = <INFILE> || die; 
		
		chomp $seq;
		chomp $quals;
		$seq = &Debarcer::trimSeq($seq, $quals);
		
		my ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq) = &Debarcer::identifyAmplicon($seq);

		next unless ( $AmpliconID eq "UNKNOWN" );
		
		if ( $AmpliconID eq "UNKNOWN" ) {
		
			my $primerList;
			($AmpliconID, $barcode_pos, $barcode, $primerList) = &Debarcer::identifyOffTarget($seq);
			$ampliconSeq = $primerList;
		
		}
		
		
		# Ignore amplicon if it is not of the expected size
		# A hack to avoid indels for now
		# next unless ( length($ampliconSeq) == $Debarcer::ampLen{$AmpliconID} );
		
		# next unless ( $barcode_pos == 0 ); # This discards about 5% of the data but barcodes not at the beginning are suspect

		print join("\t", ("OffTargetAmplicon:", $AmpliconID, $barcode_pos, $barcode, $ampliconSeq, $seq)). "\n";
		
	}

} 

close INFILE;


exit;

#########################

