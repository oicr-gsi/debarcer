

=pod

Summarize the stream from a command like

  time perl classifyReadsToAmplicons_dump.pl -fastqgz=../fastqs/Sample_01.R1.fastq.gz | 
	perl reportCommonPositionBases.pl --raw > Sample01.rawPositionalComposition.txt
  
  time cat Sample01.rawConsSeq.txt | 
	perl reportCommonPositionBases.pl --consensus > Sample01.consensusPositionalComposition.txt
  	
	
=head2 Author

Paul M Krzyzanowski
pmkrzyzanowski@gmail.com
(c) 2014-2015

=cut

# Setup general modules
# use Cwd 'abs_path';
# use File::Basename;
use Getopt::Long;

my %args = ();

GetOptions(
	"raw" => \$args{"raw"},
	"consensus" => \$args{"consensus"},
	"n=s" => \$args{"nReadsInBarcode"}
);


my %Amps = ();
my %AmpCount = ();
my $AllAmpCount = 0;
my $barcodeCutoff = 3;
if ( $args{"nReadsInBarcode"} ) {
	$barcodeCutoff = $args{"nReadsInBarcode"};
	}

print STDERR "--- Starting reportCommonPositionBases.pl ---\n";
print STDERR "reportCommonPositionBases.pl running in Raw Mode\n" if ( $args{"raw"} );
print STDERR "Using barcode depth cutoff = $barcodeCutoff\n";

while (<>) {

	chomp;

	my ($amplicon, $pos, $barcode, $seq, $barcode_count, $seqMarkup);
	
	if ( $args{"raw"} ) {
		($amplicon, $pos, $barcode, $seq) = split /\t/;
	} elsif ( $args{"consensus"} ) {
		($amplicon, $barcode, $barcode_count, $seq, $seqMarkup) = split /\t/;
		# next if ( $barcode_count < $barcodeCutoff );
		if ( $barcode_count < $barcodeCutoff ) {
			# Following touches the amplicon variables if they don't exist
			# It's useful when the depth of coverage is low and an amplicon doesn't have 
			# a deep family member to trigger indexing.
			unless ( exists $AmpCount{$amplicon} ) { 
				$AmpCount{$amplicon}->{"count"} = 0;
				my @s = split(//, $seq);
				for (my $i = 0; $i < scalar(@s); $i++) { 
					$Amps{$amplicon}->{$i} = {};
				}
			}
			next;
		}
		next if ( $amplicon eq "UNKNOWN" );
	}
	
	$AmpCount{$amplicon}->{"count"}++; # Count the instances of each amplicon exceeding the barcodeCutoff
	$AmpCount{$amplicon}->{"coverage"} += $barcode_count; # Count the read contribution of each family found which exceeds the barcodeCutoff
	$AllAmpCount++;
	
	my @s = split(//, $seq);
	for (my $i = 0; $i < scalar(@s); $i++) { 
		$Amps{$amplicon}->{$i}{$s[$i]}++;
	}
	
}


foreach my $amplicon ( sort keys %Amps ) {
	
	if ( $args{"consensus"} ) {
		printf STDERR ("%s\tdepth|count|coverage\t%d\t%d\t%d\n", $amplicon, $barcodeCutoff, $AmpCount{$amplicon}->{"count"}, $AmpCount{$amplicon}->{"coverage"});
		}
		
	foreach my $i ( sort {$a <=> $b} keys %{$Amps{$amplicon}} ) {
		
		# This is output to STDOUT and written to the *consensusPositionalComposition.bc?.txt files
		# This data is then used to generate the Amplicon Error Plots.
		# DO NOT REDIRECT TO STDERR!!
		printf ("%s\t%s", $amplicon, $i);
		my $totalReads = 0;
		foreach my $base ( qw/A C T G/ ) {
			printf("\t%s", $Amps{$amplicon}->{$i}{$base} );
			$totalReads += $Amps{$amplicon}->{$i}{$base};
			}
		printf ("\t%s", $totalReads);
		print "\n";			
			
	}
	
}

print STDERR "Total Amplicons (Including UNKNOWN if in raw mode) $AllAmpCount\n";

print STDERR "--- reportCommonPositionBases.pl Complete ---\n";
exit;
