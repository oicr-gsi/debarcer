#!/usr/bin/perl
use strict;

=pod

=head1 generateConsensusFromBAM.pl

Description

=head2 Usage

   perl generateConsensusFromBAM.pl --bam=
   
   or
   
   # perl generateConsensusFromBAM.pl --bam= --ampTable=$AMPLICON_TABLE 
   	   
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
use Data::Dumper;
use Bio::DB::Sam;  # This needs to load on a compute node.
use JSON::XS qw(encode_json decode_json);

print STDERR "--- Starting generateConsensusFromBAM.pl ---\n";
my %args = ();

GetOptions(
	"bam=s" => \$args{"bam"},
	"sampleID=s" => \$args{"sampleID"},
	"consDepth=s" => \$args{"consDepth"},  # Amplicon table filename
	"plexity=s" => \$args{"plexity"},  # Amplicon table filename
	# "strictCons" => \$args{"strictCons"},  # Strict consensus
	"downsample=s" => \$args{"downsample"},
	"justUIDdepths" => \$args{"justUIDdepths"},
	"justTargets" => \$args{"justTargets"},
	"test" => \$args{"test"}   # test mode
);

$args{"justTargets"} = 1;
my $consensusDepth = ( $args{"consDepth"} ) ? $args{"consDepth"} : 3;  # Minimum depth of a family to create a consensus call
my $nSites = ( $args{"plexity"} ) ? $args{"plexity"} : 5;  # Proxy for plexity
print STDERR "Using Consensus Depth = $consensusDepth and plexity = $nSites\n";

my $uidDepthFile = $args{"sampleID"} . ".UIDdepths.txt.gz";
open UIDDEPTHFILE, "| gzip -c > $uidDepthFile";
my $ConsensusFile = $args{"sampleID"} . ".consensusSequences.cons$consensusDepth.txt.gz";
open CONSENSUSFILE, "| gzip -c > $ConsensusFile";

# %Debarcer::primerSets = &Debarcer::listPrimers(\@ampliconsFile);
# %Debarcer::ampSeq = &Debarcer::ampliconSequences(\@ampliconsFile);
# %Debarcer::ampLen = &Debarcer::ampliconLengths(\@ampliconsFile);

my %readData = ();
my %SNVdataMaster = map { $_ => undef } qw/A C T G D I N/;  # This is a running total of all SNV types
my %familyData = ();

my $inputSeqCount = 0;
my $familySitesSeqCount = 0;
my $refgenome="/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.6.2/hg19_random.fa";
my $infile = $args{"bam"};

my %familySites = &identifyFamilySites($infile, $nSites);
my %siteAliasTable = &Debarcer::getPositionAliases("$FindBin::Bin/amplicon_tables/all_amplicons.txt");
my %ampliconInfo = ();
&Debarcer::loadAmpliconData("$FindBin::Bin/amplicon_tables/all_amplicons.txt", \%ampliconInfo);
my %invalidBarcodes = &loadBarcodeMaskFile($args{"sampleID"});
# print Dumper(\%ampliconInfo);
# print Dumper(\%siteAliasTable);
# print Dumper(\%invalidBarcodes);

my $sam = Bio::DB::Sam->new(-bam => $infile, -fasta => $refgenome );
my $bam = $sam->bam();
my $header = $bam->header();

my $iterator = $sam->features(-iterator=>1);
while (my $alignment = $iterator->next_seq) { 

	next if ( $args{"test"} & (rand() < 0.001) );  # For sampling.
	
	$inputSeqCount++;
	
	my $chrom = $alignment->seq_id;
	my $chromStart = $alignment->start();
	
	# Check that this read is is the set of analyzable amplicons
	my $AmpliconID = "$chrom:$chromStart";
	next unless ( exists $familySites{$AmpliconID} );
	$familySitesSeqCount++;
	
	my $read_dna = $alignment->query->dna();
	my ($bc_position, $barcode) = &Debarcer::extractBarcodeQuick($read_dna);
	
	next unless ( $bc_position == 0 );

	# Skip the barcode if it's in the mask hash loaded from the maskfile
	if ( exists $invalidBarcodes{$AmpliconID}->{$barcode} ) {
		# print "Skipping invalidBarcodes{$AmpliconID}->{$barcode}\n"; 
		next unless ( $args{"justUIDdepths"} ); # Do not skip an invalid barcode if we're only generating the UID depth file.
	}
	
	# Count this instance of the barcode family
	$familyData{$AmpliconID}->{$barcode}{"count"}++;
	next if ( $args{"justUIDdepths"} );
	
	# Determine the base call or insertion/deletion status wrt to the start of the alignment, i.e.
	# chromosomal position given in $chromStart.  Therefore $basecalls[0] is for $chromStart + 0
	my @basecalls = &calculateBasecalls($alignment);

	# Save the raw base calls for future consensus calling
	# Write the base by base calls to a file....
	push(@{$familyData{$AmpliconID}->{$barcode}{"raw"}}, join("", @basecalls));

	# Store the individual basecalls in the readData hash
	for ( my $i = 0; $i < scalar(@basecalls); $i++ ) {
		$basecalls[$i] =~ tr/acgt/ACGT/;
		$readData{$AmpliconID}->{$i}{$basecalls[$i]}++;
		$SNVdataMaster{$basecalls[$i]}++;
	}
	

}

# print Dumper(\%SNVdataMaster);

# Just print the UID depth file if the --justUIDdepths flag is set.
if ( $args{"justUIDdepths"} ) {
	foreach my $amp ( keys %familyData ) {
		foreach my $barcode ( keys %{$familyData{$amp}} ) {
			( my $printableAmpID = $amp ) =~ s/\t/\:/;
			printf UIDDEPTHFILE ("%s\t%s\t%s\n", $printableAmpID, $barcode, $familyData{$amp}->{$barcode}{"count"} );
		}
	}
	exit;
}

# Compute consensus sequences
my %consReadData = ();

foreach my $amp ( keys %familyData ) {
	my ($AmpliconCount, $AmpliconCoverage) = (0, 0);
	( my $printableAmpID = $amp ) =~ s/\t/\:/;
	foreach my $barcode ( keys %{$familyData{$amp}} ) {
		
		if ( $familyData{$amp}->{$barcode}{"count"} >= $consensusDepth ) {
			$AmpliconCount++;
			$AmpliconCoverage += $familyData{$amp}->{$barcode}{"count"};
		}
		
		$familyData{$amp}->{$barcode}{"consensus"} = &generateConsensus(@{$familyData{$amp}->{$barcode}{"raw"}}, $consensusDepth);
		# $familyData{$amp}->{$barcode}{"consensus"} = &generatePhyloConsensus(@{$familyData{$amp}->{$barcode}{"raw"}}, $consensusDepth);  # Still in progress.

		# Write the UID depth information to a file
		printf UIDDEPTHFILE ("%s\t%s\t%s\n", $printableAmpID, $barcode, $familyData{$amp}->{$barcode}{"count"});
		printf CONSENSUSFILE ("%s\t%s\t%s\t%s\n", $printableAmpID, $barcode, $familyData{$amp}->{$barcode}{"count"}, $familyData{$amp}->{$barcode}{"consensus"} );
		
		my @cons = split(//, $familyData{$amp}->{$barcode}{"consensus"});
		for (my $i = 0; $i < scalar(@cons); $i++) {
			$cons[$i] =~ tr/acgt/ACGT/;
			$consReadData{$amp}->{$i}{$cons[$i]}++;
		}
	}
	print STDERR "$printableAmpID\tdepth|count|coverage\t$consensusDepth\t$AmpliconCount\t$AmpliconCoverage\n";
}

# Consider storing family Data at this point
# See http://www.perlmonks.org/?node_id=510202
# my $JSONdata = encode_json(\%familyData);
# open(OFIL,">fileHash.dat"); print OFIL $JSONdata; close(OFIL);
# print Dumper(\%familyData);


# output some data
my @SNVtypes = sort keys %SNVdataMaster;
print join ("\t", "#AmpliconChromStart", "Alias", "Position", "ProbableRef", 
	"raw" . join("\traw", @SNVtypes, "Depth"), 
	"cons" . join("\tcons", @SNVtypes, "Depth"), 
	"\n");

# Identify a minimal depth to report positions within amplicons. 
# This eliminates the reporting of bases indexed from secondary PCR products
my %ampliconReportingThreshold = &calculateAmpliconReportingThreshold(\%readData);
# print Dumper(\%ampliconReportingThreshold);

# print Dumper(\%ampliconInfo);  # 	
	
foreach my $amp ( keys %readData ) {
	# print Dumper(\%{$readData{$amp}});
	foreach my $position ( sort { $a <=> $b } keys %{$readData{$amp}} ) {

		# output raw data
		my ($probableBase, $n, $rawDataString) = ('', 0, '');
		my $d = 0;  #Depth
		foreach ( @SNVtypes ) {
			# Save the base identity if the count is higher than any previous observed count
			if ( $readData{$amp}{$position}{$_} > $n ) {
				$probableBase = $_;
				$n = $readData{$amp}{$position}{$_};
			}
			$rawDataString .= "\t" . $readData{$amp}{$position}{$_};
			$d += $readData{$amp}{$position}{$_};
		}

		# Do not report this site if it's a low coverage position
		next if ( $d < $ampliconReportingThreshold{$amp} );
		
		my $ampliconName = $siteAliasTable{$amp};
		my ($chrom, $chromStart) = split(/:/, $amp);
		my $genomePosition = $chromStart + $position;

		# This section will restrict reporting of amplicon positions to those within
		# the target window specified in the all_amplicons.txt file
		if ( $args{"justTargets"} ) {
			# If a target window exists, skip $genomePosition unless it falls within start and end
			if ( exists $ampliconInfo{$ampliconName}{"TargetWindow"} ) {
				unless ( $genomePosition >= $ampliconInfo{$ampliconName}{"TargetWindow"}{"start"} & $genomePosition <= $ampliconInfo{$ampliconName}{"TargetWindow"}{"end"} ) {
					# print STDERR "Skipping $ampliconName $amp Position $position $genomePosition\n";
					next;
				}
			}
			# To gen here, no target window exists, so report everything
		}


		printf("%s\t%s\t%s",
			$amp,
			$ampliconName,
			$position);
		print "\t$probableBase";
		print $rawDataString;
		print "\t$d";

		# output consensus data
		$d = 0;
		foreach ( @SNVtypes ) {
			print "\t" . $consReadData{$amp}{$position}{$_};
			$d += $consReadData{$amp}{$position}{$_};
		}
		print "\t$d";

		print "\n";

	}
}


print STDERR "Raw reads read from $infile: $inputSeqCount\n";
print STDERR "Raw reads in family sites read from $infile: $familySitesSeqCount\n";

print STDERR "--- generateConsensusFromBAM.pl Complete ---\n";

close UIDDEPTHFILE;
close CONSENSUSFILE;

exit;

#################################################################################################


sub identifyFamilySites {
	my $inBam = shift @_;
	my $nSites = shift @_;
	my %familySites = ();
	
	my @allSites = `samtools view -s 0.1 $inBam | cut -f 3,4`;
	foreach my $site ( @allSites ) { chomp $site; $site =~ s/\t/:/; $familySites{$site}++; }
	my @goodSites = (sort { $familySites{$b} <=> $familySites{$a} } keys %familySites);
	@goodSites = @goodSites[0..$nSites];  # Take the top n-1
	print STDERR "Compiling info for Family Sites:\n";
	for ( my $site = 0; $site < scalar( @goodSites ); $site++ ) {
		print STDERR "$goodSites[$site]\t" . $familySites{$goodSites[$site]} . "\n";
		splice( @goodSites, $site, 1) if ( $site =~ /\*/ );  #Remove the current site if it's unmapped.
		}
	print STDERR "Note:  If present, the '* 0' (i.e. unmapped) site has been dropped\n";
	
	my %goodHash = ();
	@goodHash{@goodSites} = 1;
	return %goodHash;
	
}

sub calculateBasecalls {

=pod

=head1 calculateBase calls

A function that returns the position by position calls [ACTGDI] reported by a
Bio::DB::Bam::Alignment object

=cut

	my $alignment = shift @_;
	my @basecalls = '';
	my $debug = 0;
	
	my $chromStart = $alignment->start;
	
	# For testing
	if ( $debug ) {
		my $CIGAR = $alignment->cigar_str;
		return 1 unless ( $CIGAR =~ /[ID]/ );
		print "$CIGAR\n";
	}
	
	my ($ref,$matches,$query) = $alignment->padded_alignment;
	# print "Ref:  $ref\n      $matches\nRead: $query\n\n";
	
	# Since query is longer than ref due to adapters, etc.,
	# trim the sequences
	$ref =~ /^(-+).+?(-+)$/;
	
	($ref, $query) = ( substr($ref, length($1), length($ref) - (length($1 . $2)) ), substr($query, length($1), length($query) - (length($1 . $2)) ) );
	($matches) = ( substr($matches, length($1), length($matches) - (length($1 . $2)) ) );
	print "Mismatch:\n" if ( $debug & ($ref ne $query) );
	print "Ref:  $ref\n      $matches\nRead: $query\n" if ( $debug );
	
	# <STDIN>;  # Pause for keypress
	
	# Index the sequence by genomic position
	my $genomicOffset = 0;
	my $inInsertion = 0;
	
	my @Ar = split(//, $ref);
	my @Aq = split(//, $query);
		
	for (my $i = 0; $i < scalar(@Ar); $i++) {
		if ( $Ar[$i] eq $Aq[$i] ) {
			$basecalls[$genomicOffset] = $Ar[$i];  # This is a match
			$genomicOffset++;
			$inInsertion = 0 if ( $inInsertion );
		} elsif ( $Aq[$i] eq '-' ) {
			$basecalls[$genomicOffset] = 'D';  # This is deletion of a ref base
			$genomicOffset++;
			$inInsertion = 0 if ( $inInsertion );
		} elsif ( $Ar[$i] eq '-' ) {
			next if ( $inInsertion );
			$basecalls[($genomicOffset-1)] = 'I';  # This is an insertion at the genomic position prior to the first position of the insertion
			$inInsertion = 1;
		} elsif ( $Ar[$i] ne $Aq[$i] ) {
			$basecalls[$genomicOffset] = $Aq[$i];  # This is just a mismatch
			$basecalls[$genomicOffset] =~ tr/ACTG/actg/;  # Convert to lowercase to symbolize a non-reference base
			$genomicOffset++;
			$inInsertion = 0 if ( $inInsertion );
		} else {
			print join(" ", $Ar[$i], $Aq[$i], $i, "There is a problem.", "\n");
		}
			
	}
	
	print join("", "CIGAR:", @basecalls, "*\n\n") if ( $debug );
	
	return @basecalls;
			
}

sub generateConsensus {

=pod

Function to return a consensus sequence given an array of CIGAR-like alignment strings.

=cut
	
	my $minDepth = pop @_;
	my @rawSeqs = @_;
	return if ( scalar(@rawSeqs) < $minDepth );  # Return a blank consensus if one can't be called given minDepth

	my $cons;
	my %rawBases = ();

	# index by position in the AoA
	foreach my $seq ( @rawSeqs ) {
		my @s = split(//, $seq);
		for ( my $i = 0; $i < scalar(@s); $i++) {
			$rawBases{$i}{$s[$i]}++;
		}
	}

	# print Dumper(\%rawBases);
	
	foreach my $i ( sort {$a <=> $b} keys %rawBases ) {
		my @basesHere = sort { $rawBases{$i}{$b} <=> $rawBases{$i}{$a} } keys %{$rawBases{$i}};  # sort bases by descending count
		my $commonBase = $basesHere[0];
		
		# if commonBase is in [acgtDI] it must be highly abundant.  There should be a test here.
		if ( $commonBase =~ /[acgtDI]/ ) {
			# print "## Checking this non-reference base: $commonBase :" . Dumper(\%{$rawBases{$i}});
			my $nCommonBase = $rawBases{$i}{$commonBase};
			my $depthHere = 0;
			$depthHere += $rawBases{$i}{$_} foreach ( @basesHere );
			if ( $depthHere <= 20 ) {
				unless ( $nCommonBase == $depthHere ) {
					# print "Changing $commonBase to $basesHere[1] because of non-unanimous evidence: " . Dumper(\%{$rawBases{$i}});
					$commonBase = $basesHere[1];
				}
			} else {
				my $alleleRatio = ($nCommonBase / $depthHere) ;
				unless ( $alleleRatio >= 0.90 ) {
					# print "Changing $commonBase to $basesHere[1] because of $alleleRatio: " . Dumper(\%{$rawBases{$i}});
					$commonBase = $basesHere[1];
				}
			}
		}
		
		$cons .= $commonBase;
	}

	return $cons;
}

sub generatePhyloConsensus {

=pod

Function to return a consensus sequence, based on phylogenetic relationships, given an array of CIGAR-like alignment strings.

=cut
	
	my $minDepth = pop @_;
	my @rawSeqs = @_;
	return if ( scalar(@rawSeqs) < $minDepth );  # Return a blank consensus if one can't be called given minDepth

	my $cons;
	
	my %seqHash = ();
	$seqHash{$_}++ foreach ( @rawSeqs );
	
	print Dumper(\%seqHash);

	
	
	
	return $cons;
}

sub loadBarcodeMaskFile {

=pod

Load the $sample.barcode_mask file, if it exists

=cut

	my $sampleID = shift @_;
	my %barcodeMask = ();
	my $maskfile = "$sampleID.barcode_mask";
	
	if ( -e $maskfile ) {
		# print STDERR "Loading $maskfile\n";
		open INFILE, $maskfile;
		while (<INFILE>) {
			next unless (/Mask/);
			my @line = split(/\t/);
			shift @line;
			my $amp = shift @line;
			my $invalidBarcode = shift @line;
			$amp =~ s/\:/\t/;
			$barcodeMask{$amp}->{$invalidBarcode}++;
			}
		close INFILE;
	}

	return %barcodeMask;
}

sub calculateAmpliconReportingThreshold {

=pod

=head2 calculateAmpliconReportingThreshold

This function identifies a minimal threshold for each amplicon,
below which positions aren't reported in the cons<depth>.txt files

=cut

	my $href = shift @_;  # This is %readData
	# readData format is
	# $readData{$AmpliconID}->{$position}{$basecall};
	my %depthCuts = ();
	
	# find the maximum depth for each amplicon ( usually the first site )
	foreach my $amplicon ( keys %$href ) {
		foreach my $position ( keys %{$href->{$amplicon}} ) {
			my $depthHere = 0;
			$depthHere += $href->{$amplicon}{$position}{$_} foreach ( keys %{$href->{$amplicon}{$position}} );
			$depthCuts{$amplicon} = $depthHere if ( $depthHere > $depthCuts{$amplicon} );
		}
	}
	
	# Adjust depth cuts downward to a percentage of max
	foreach my $amplicon ( keys %depthCuts ) { $depthCuts{$amplicon} = int( $depthCuts{$amplicon} * 0.1 ); }
	
	return %depthCuts;
}

# SVNID: $Id: generateConsensusFromBAM.pl 387 2016-02-09 22:08:35Z pkrzyzanowski $
