#!/usr/bin/perl
use strict;

=pod

=head1 classifyReadsToAmplicons.pl

Description

=head2 Usage

   perl classifyReadsToAmplicons.pl --fastqgz=../fastqs/Sample_01.R1.fastq.gz --ampTable=$AMPLICON_TABLE 
   	   
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
use lib '/u/pkrzyzanowski/perl/usr/local/lib/perl/5.10.1/';
# use Bio::Grep;
# use String::Approx 'asubstitute';
# use re::engine::PCRE;

print STDERR "--- Starting classifyReadsToAmplicons.pl ---\n";
my %args = ();

GetOptions(
	"fastqgz=s" => \$args{"fastqgz"},
	"bamfile=s" => \$args{"bamfile"},
	"plexity=s" => \$args{"plexity"},      # The number of amplicons expected
	"ampTable=s" => \$args{"ampTable"},  # Amplicon table filename
	"strictCons" => \$args{"strictCons"},  # Strict consensus
	"dump" => \$args{"dump"},  # dump mode
	"consensus" => \$args{"consensus"},  # consensus mode
	"test" => \$args{"test"}   # test mode
);

# Put the amplicon table into a global array
unless ( $args{"ampTable"} ) { # If no argument given use this default file
	$args{"ampTable"} =  "$ENV{BHOME}/amplicon_tables/all_amplicons.txt"; 
	}
open AMPS, $args{"ampTable"} or die "Can't open table of amplicon information\n";
my @ampliconsFile = <AMPS>;
# print STDERR "Amplicons:\n@ampliconsFile\n";
close AMPS;

$args{"plexity"} = 5; # Default is to expect 5 amplicons.

%Debarcer::primerSets = &Debarcer::listPrimers(\@ampliconsFile);
%Debarcer::ampSeq = &Debarcer::ampliconSequences(\@ampliconsFile);
%Debarcer::ampLen = &Debarcer::ampliconLengths(\@ampliconsFile);

my %readData = ();
my $inputSeqCount = 0;


if ( $args{"dump"} && $args{"fastqgz"} ) { # Run dump mode

	my $infile = $args{"fastqgz"};
	if ( $args{"test"} ) {
		open INFILE, "gunzip -c $infile | head -n 400000 |";
		} else {
		open INFILE, "gunzip -c $infile |";
		}

	while (my $header = <INFILE>) {
		my $seq = <INFILE> || die; 
		next if ($header =~ /\+/); # Skip since $seq is the string of Q-values.
		my $separator = <INFILE> || die; # Read the + symbol
		my $quals = <INFILE> || die; 
		$inputSeqCount++;
		
		chomp $seq;
		chomp $quals;
		$seq = &Debarcer::trimSeq($seq, $quals);
		
		my ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq) = &Debarcer::identifyAmplicon($seq);

		# Ignore amplicon if it is not of the expected size
		# A hack to avoid indels for now
		unless ( length($ampliconSeq) == $Debarcer::ampLen{$AmpliconID} ) {
#			print STDERR join("\t", ("DiscardLengthMismatch:", $AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
			next;
			}
		
		# This discards about 5% of the data but barcodes not at the beginning are suspect
		unless ( $barcode_pos == 0 ) {
#			print STDERR join("\t", ("DiscardBarcodePosition:", $AmpliconID, $barcode_pos, $barcode, $seq)). "\n";
			next;
			}

		print join("\t", ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
		
	}

	print STDERR "Raw reads read from $infile: $inputSeqCount\n";
	close INFILE;

} elsif ( $args{"dump"} && $args{"bamfile"} ) { # Run dump mode

	# This hashes keys contain the most abundant start sites according to the plexity
	my %validAmplicons = &Debarcer::detectAmplicons($args{"bamfile"}, $args{"plexity"});

	my $infile = $args{"bamfile"};
	if ( $args{"test"} ) {
		open INFILE, "samtools view -s 0.001 -X $infile |"; # Sample a fraction of the file.
		} else {
		open INFILE, "samtools view -X $infile |";
		}

	while (<INFILE>) {
		
		my @line = split /\t/;
		my $flag = $line[1];
		my $chrom = $line[2];
		next if ( $chrom eq "*" );
		my $chromStart = $line[3];
		my $CIGAR = $line[5];
		my $seq = $line[9];
		$inputSeqCount++;
		
		# For a bam file, amplicon ID becomes the start site of the alignment.
		my $_AmpliconID = "$chrom:$chromStart";
		next unless exists ( $validAmplicons{$_AmpliconID} );
		
		chomp $seq;
		# not needed for bam files
		# chomp $quals;
		# $seq = &Debarcer::trimSeq($seq, $quals);
		
		my ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq) = &Debarcer::identifyAmpliconBam($seq, $CIGAR);
		

		
		# Ignore amplicon if it is not of the expected size
		# A hack to avoid indels for now
		unless ( length($ampliconSeq) == $Debarcer::ampLen{$AmpliconID} ) {
#			print STDERR join("\t", ("DiscardLengthMismatch:", $AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
			# next;  
			}
		
		# This discards about 5% of the data but barcodes not at the beginning are suspect
		unless ( $barcode_pos == 0 ) {
#			print STDERR join("\t", ("DiscardBarcodePosition:", $AmpliconID, $barcode_pos, $barcode, $seq)). "\n";
			next;
			}

		print join("\t", ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
		
	}
	
	print STDERR "Raw reads read from $infile: $inputSeqCount\n";
	close INFILE;

} elsif ( $args{"consensus"} && $args{"fastqgz"}  ) {
	
	my $infile = $args{"fastqgz"};
	if ( $args{"test"} ) {
		open INFILE, "gunzip -c $infile | head -n 400000 |";
		} else {
		open INFILE, "gunzip -c $infile |";
		}

	while (my $header = <INFILE>) {
		my $seq = <INFILE> || die; 
		next if ($header =~ /\+/); # Skip since $seq is the string of Q-values.
		chomp $seq;
		$inputSeqCount++;
		
		my ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq) = &Debarcer::identifyAmplicon($seq);

		# Ignore amplicon if it is not of the expected size
		# A hack to avoid indels for now
		unless ( length($ampliconSeq) == $Debarcer::ampLen{$AmpliconID} ) {
#		print STDERR join("\t", ("Amplicon_Discarded:", $AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
			next;
			}
		
		# This discards about 5% of the data but barcodes not at the beginning are suspect
		unless ( $barcode_pos == 0 ) {
#			print STDERR join("\t", ("Amplicon_Discarded:", $AmpliconID, $barcode_pos, $barcode, $seq)). "\n";
			next;
			}

		# print join("\t", ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
		

		# Store the data in the readData hash
		$readData{$AmpliconID}->{$barcode}{"count"}++;
		my @bases = split(//, $ampliconSeq);
		for ( my $i = 0; $i < scalar(@bases); $i++ ) {
			$readData{$AmpliconID}->{$barcode}{"PSWM"}{$i}{$bases[$i]}++;	
		}
		
	}

	print STDERR "Raw reads read from $infile: $inputSeqCount\n";
	close INFILE;

} elsif ( $args{"consensus"} && $args{"bamfile"}  ) {

	# This hashes keys contain the most abundant start sites according to the plexity
	my %validAmplicons = &Debarcer::detectAmplicons($args{"bamfile"}, $args{"plexity"});
	
	my $infile = $args{"bamfile"};
	if ( $args{"test"} ) {
		open INFILE, "samtools view -s 0.01 -X $infile |"; # Sample a fraction of the file.
		} else {
		open INFILE, "samtools view -X $infile |";
		}

	while (<INFILE>) {
		
		my @line = split /\t/;
		my $flag = $line[1];
		my $chrom = $line[2];
		next if ( $chrom eq "*" );
		my $chromStart = $line[3];
		my $CIGAR = $line[5];
		my $seq = $line[9];
		$inputSeqCount++;
		
		# For a bam file, amplicon ID becomes the start site of the alignment.
		my $_AmpliconID = "$chrom:$chromStart";
		next unless exists ( $validAmplicons{$_AmpliconID} );
				
		chomp $seq;
		# not needed for bam files
		# chomp $quals;
		# $seq = &Debarcer::trimSeq($seq, $quals);
		
		my ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq) = &Debarcer::identifyAmpliconBam($seq, $CIGAR);
		
		# For a bam file, amplicon ID becomes the start site of the alignment.
		$AmpliconID = "$chrom:$chromStart";
		
		# Ignore amplicon if it is not of the expected size
		# A hack to avoid indels for now
		# unless ( length($ampliconSeq) == $Debarcer::ampLen{$AmpliconID} ) {
#			print STDERR join("\t", ("DiscardLengthMismatch:", $AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
			# next;  
			# }
		
		# This discards about 5% of the data but barcodes not at the beginning are suspect
		next unless ( $barcode );

		# print join("\t", ($AmpliconID, $barcode_pos, $barcode, $ampliconSeq)). "\n";
		
		# Store the data in the readData hash
		$readData{$AmpliconID}->{$barcode}{"count"}++;
		my @bases = split(//, $ampliconSeq);
		my $nBases = scalar(@bases);
		
		# The first sequence sets the "real" number of bases.
		# This is a hack to avoid indel containing sequences
		next if ( exists $readData{$AmpliconID}->{$barcode}{"PSWM"} && ( $nBases != scalar keys %{$readData{$AmpliconID}->{$barcode}{"PSWM"}} ) );
		
		for ( my $i = 0; $i < scalar(@bases); $i++ ) {
			$readData{$AmpliconID}->{$barcode}{"PSWM"}{$i}{$bases[$i]}++;	
		}
			
	}
	
	print STDERR "Raw reads read from $infile: $inputSeqCount\n";
	close INFILE;

}

# output some data
foreach my $amp ( keys %readData ) {
	foreach my $bc ( keys %{$readData{$amp}} ) {
		my $n = $readData{$amp}{$bc}{"count"};
		
		my $consensusSeqLong = &Debarcer::computeConsensusLong($readData{$amp}->{$bc}{"PSWM"});
		
		my $consensusSeq;
		if ( $args{"strictCons"} ) {
			$consensusSeq = &Debarcer::computeConsensus_strict($readData{$amp}->{$bc}{"PSWM"}, $amp);
			} else {
			$consensusSeq = &Debarcer::computeConsensus($readData{$amp}->{$bc}{"PSWM"});
			}
		
		# printf("%s\t%s\t%s\t%s\t%s\n", $amp, $bc, $n, $consensusSeq, $consensusSeqLong) if ( $n > 1 );
		printf("%s\t%s\t%s\t%s\t%s\n", $amp, $bc, $n, $consensusSeq, $consensusSeqLong);
		
	}
}


print STDERR "--- classifyReadsToAmplicons.pl Complete ---\n";
exit;

#########################
