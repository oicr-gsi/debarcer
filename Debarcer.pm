#!/usr/bin/perl
package Debarcer;

use strict;

our %ampSeq = ();
our %ampLen = ();
our %primerSets = ();

use lib "$ENV{'BHOME'}/utils/Text-Levenshtein-0.12/lib/";
use Text::Levenshtein qw(distance);
use Data::Dumper; # for logging calls
use Log::Log4perl; # see http://search.cpan.org/~mschilli/Log-Log4perl-1.46/lib/Log/Log4perl.pm

# FIXME: is this where the logging config file should go?
Log::Log4perl::init_and_watch("$ENV{'BHOME'}/config/log4perl.conf");

# $BHOME/config/log4perl.conf looks something like this:
#--------------------
# log4perl.logger.debarcer            = DEBUG, FileAppender
# log4perl.appender.FileAppender      = Log::Log4perl::Appender::File
# log4perl.appender.FileAppender.filename = debarcer.log
# log4perl.appender.FileAppender.layout   = Log::Log4perl::Layout::SimpleLayout
#--------------------

sub identifyAmpliconBam {
	my $seq = shift @_;
	my $CIGAR = shift @_;
	my $whichAmplicon = 'UNKNOWN';
	my $ampliconSeq = '';

	my $fuzzyMatching = 1;
	# my $fuzzyMatching = 0;
	
	my ($bc_position, $barcode, $isRev) = &extractBarcodeQuick($seq);
	
	# Report unknown if there was no valid barcode
	return ($whichAmplicon, $bc_position, $barcode, $ampliconSeq) unless ( $bc_position == 0 );	
	
	$CIGAR =~ m/^(\d+)S.+(\d+)S$/g;
	$ampliconSeq = substr($seq, $1, length($seq) - $1 - $2);
	$ampliconSeq = &revcom($ampliconSeq) if ( $isRev );
	$whichAmplicon = "";  # The amplicon name will become the genomic start of the CIGAR alignment
	
	return ($whichAmplicon, $bc_position, $barcode, $ampliconSeq);	
} # End sub identifyAmpliconBAM


sub identifyAmplicon {
	my $seq = shift @_;

	my $fuzzyMatching = 1;
	# my $fuzzyMatching = 0;

	my $logger = Log::Log4perl->get_logger('debarcer');
	
	my ($bc_position, $barcode) = &extractBarcode($seq);
	my $whichAmplicon = 'UNKNOWN';
	my $ampliconSeq = '';
	
	foreach my $ampID ( keys %primerSets ) {
		my ($fwd, $rev) = @{$primerSets{$ampID}};
		
		if ($fuzzyMatching) {

			my $fwdMatch = &fuzzyMatch2($fwd, $seq);
			next unless ( $fwdMatch );
			my $revMatch = &fuzzyMatch2($rev, $seq);
			next unless ( $revMatch );
			
			$seq =~ m/$fwdMatch([ACTG]+?)$revMatch/g;
			$ampliconSeq = $1;
			$whichAmplicon = $ampID;
			last;
			
		} else {

			if ( $seq =~ $fwd && $seq =~ $rev ) {
				$seq =~ m/$fwd([ACTG]+?)$rev/g;
				$ampliconSeq = $1;
				$whichAmplicon = $ampID;
				last;
			}		

		}
			
	}

	$logger->debug("\n$seq\nAMPLICON= $whichAmplicon $barcode $ampliconSeq\n");
	
	return ($whichAmplicon, $bc_position, $barcode, $ampliconSeq);	
}
=pod

=head fuzzyMatch

Return a string that matches between a target and query,
allowing up to one mismatch

=cut

sub fuzzyMatch {
	my ($query, $subject) = @_;
	my $logger = Log::Log4perl->get_logger('debarcer');

	return $query if ( $subject =~ /$query/ ); # If an exact match exists return the original query

	# Test for a partial match to the current sequence prior to doing the expensive regex search
	my $halfA = substr($query, 0, 10);
	my $halfB = substr($query, 10);
	if ( $subject =~ /$halfA/ || $subject =~ /$halfB/ ) {
	
		for ( my $i = 0; $i < length($query); $i++ ) {
			my $newQuery = $query;
			substr($newQuery, $i, 1) = '[ACGT]?'; # match 1 bp deletion or a substitution at position $i
			$logger->debug("$query $newQuery\n");
			return $newQuery if ( $subject =~ /$newQuery/ ); # If the regex matches return the regex in newQuery
		}
		
	}
	
	return 0;
}


=pod

=head fuzzyMatch2

Return a string that matches between a target and query,
allowing up to two mismatches

=cut

sub fuzzyMatch2 {
	my ($query, $subject) = @_;
	my $logger = Log::Log4perl->get_logger('debarcer');

	return $query if ( $subject =~ /$query/ ); # If an exact match exists return the original query

	# Test for a partial match to the current sequence prior to doing the expensive regex search
	my $halfQueryPos = int(length($query) / 2);
	my $halfA = substr($query, 0, $halfQueryPos);
	my $halfB = substr($query, $halfQueryPos);
	if ( $subject =~ /$halfA/ || $subject =~ /$halfB/ ) {
	
		for ( my $i = 0; $i < length($query); $i++ ) {
			my $newQuery = $query;
			substr($newQuery, $i, 1) = '[ACGT]?'; # match 1 bp deletion or a substitution at position $i
			$logger->debug("$query $newQuery\n");
			return $newQuery if ( $subject =~ /$newQuery/ ); # If the regex matches return the regex in newQuery
		}
		
		# Match two matches without deletion
		for ( my $i = 0; $i < (length($query) - 1); $i++ ) {
			for ( my $j = $i+1; $j < length($query); $j++ ) {
				my $newQuery = $query;
				substr($newQuery, $i, 1) = '.'; # match 1bp wildcard at position $i
				substr($newQuery, $j, 1) = '.'; # match 1bp wildcard at position $j
				$logger->debug("$query $newQuery\n");
				return $newQuery if ( $subject =~ /$newQuery/ ); # If the regex matches return the regex in newQuery
			}
		}
		
	}
	
	return 0;
}



sub computeConsensus {
	my $href = shift @_;
	my @cons = '';
	my $logger = Log::Log4perl->get_logger('debarcer');
	
	$logger->debug(Dumper($href));

	foreach my $i ( sort {$a <=> $b} keys %$href ) {
		$logger->debug("$i " . $href->{$i} . "\n");
		
		# lazy but works
		my $commonBase = (sort {$href->{$i}{$b} <=> $href->{$i}{$a}} keys %{$href->{$i}})[0];
		
		push(@cons, $commonBase);		
		
	}
	
	my $c = join("", @cons);
	
	return $c;

} # end sub computeConsensus


sub computeConsensusFromSeqarray {
	my $aref = shift @_;
	my $n_seqs = scalar @{$aref};
	my @cons = '';
	my $logger = Log::Log4perl->get_logger('debarcer');
		
	$logger->debug(Dumper($aref));
	
	# Determine most likely length of sequence
	my %CIGARs = ();
	foreach my $set ( @{$aref} ) {
		$CIGARs{$set->[2]}++;
	}
	
	if ( (scalar keys %CIGARs > 1) & 0 ) {  # <--- Set 0/1 for printing.
		printf("NSeqs=%d NumCIGARs=%s\n", $n_seqs, scalar keys %CIGARs);
		$logger->debug(Dumper($aref));
		$logger->debug(Dumper(\%CIGARs));
	}

	# Build the consensus
	#
	# Build the Position Weight Matrix
	# FIXME: break the body of this loop into a couple of subroutines to improve readability
	my %readData = ();
	foreach my $set ( @{$aref} ) {
		my ($seq, $quals, $CIGAR) = @$set;

		$logger->debug(join("\t", @$set) . "\n");
		$logger->debug(Dumper($aref));
			
		# Correct the seq for indels reported by the CIGAR string
		next if ( $CIGAR =~ /\dI/ ); # Skip insertions
		if ( $CIGAR =~ /D/ ) { # Edit sequence to include indels
			my $D_index = 0;
			while ( $CIGAR =~ m/([0-9]+)([MIDNSHPX=])/g ) {
				$D_index += $1;
				if ( $2 eq "D" ) {
					# FIXME: not sure what you want logged here
					my $test = 0;
					print "\n\nEditing using $CIGAR:\n" if ($test);
					print "$seq\n" if ($test);
					substr($seq, (-1 * $1)) = '';
					print "$seq\n" if ($test);
					substr($seq, $D_index, 0, "-" x $1);
					print "$seq\n" if ($test);
					# <STDIN>;
				}
			}
		}
		
		# Remove first 25-26 bases consisting of barcode and hairpin primer
		# This helps reduce reported errors in the long form consensus seqs.
		pos $CIGAR = 0;
		$CIGAR =~ m/([0-9]+)S/g;
		substr($seq, 0, $1, '');
		substr($quals, 0, $1, '');
				
		my @bases = split(//, $seq); 
		for ( my $i = 0; $i < scalar(@bases); $i++ ) {
			$readData{$i}->{$bases[$i]}++;	
		}
	}

	$logger->debug(Dumper(\%readData));

	# Identify common bases.
	my @cons = '';
	my @consLong = '';
	foreach my $i ( sort {$a <=> $b} keys %readData ) {
		
		# lazy but works
		my $commonBase = (sort { $readData{$i}{$b} <=> $readData{$i}{$a} } keys %{$readData{$i}})[0];
		push(@cons, $commonBase);		
		
		# Compute long form consensus
		my @bases = keys %{$readData{$i}};
		my $commonBaseLong = '';
		if ( scalar @bases > 1 ) {
			$commonBaseLong = join("", ("[", %{$readData{$i}}, "]"));
		} else {
			$commonBaseLong = $bases[0]; # Should be only one base! FIXME: put in an assertion to check this...
		}
		push(@consLong, $commonBaseLong);		
				
	}
	
	my $c = join("", @cons);
	my $cLong = join("", @consLong);
	
	return ($c, $cLong);

} # end sub computeConsensusFromSeqarray



sub computeConsensus_strict {
	my $href = shift @_;
	my $ampliconID = shift @_;
	
	my $ampliconSeq = $ampSeq{$ampliconID};
	my @cons = '';
	my $logger = Log::Log4perl->get_logger('debarcer');
	
	$logger->debug(Dumper($href));

	foreach my $i ( sort {$a <=> $b} keys %$href ) {
				
		my @bases = keys %{$href->{$i}};
		my $commonBase = '';
		if ( scalar @bases > 1 ) {

			# get the reference base from the amplicon definition
			my $refBase = substr($ampliconSeq, $i, 1);
			my $altBase = '';
			
			# calculate the total depth of reads at this position
			my $depth = 0;
			foreach ( values %{$href->{$i}} ) { $depth += $_; }

			# FIXME: convert this to logging call
			# printf ("Before Switch altBase? %s %s %s %s %s\n", $ampliconID, $i, $refBase, join("", ("[", %{$href->{$i}}, "]")), $depth);
					
			# scan through bases and assign altBase to a value if an
			# alternate allele frequency about 80% is found
			foreach my $base ( @bases ) {
				if ( $base ne $refBase ) {
					my $altBaseFreq = ( $href->{$i}{$base} / $depth );
					# FIXME: convert this to logging call
					# printf ("Switch altBase? %s %s %s %s\n", $refBase, $base, $depth, $altBaseFreq);
					$altBase = $base if ( $altBaseFreq > 0.8 );
				}
			}
			
			$commonBase = ( $altBase ) ? $altBase : $refBase;
				
		} else {
			$commonBase = $bases[0]; # Should be only one base!
		}
				
		push(@cons, $commonBase);		
		
	}
	
	my $c = join("", @cons);
	
	return $c;

} # end sub computeConsensus



sub computeConsensusLong {
	my $href = shift @_;
	my @cons = '';
	my $logger = Log::Log4perl->get_logger('debarcer');
	
	$logger->debug(Dumper($href));

	foreach my $i ( sort {$a <=> $b} keys %$href ) {
		$logger->debug("$i " . $href->{$i} . "\n");
		
		# lazy but works
		# my $commonBase = (sort {$href->{$i}{$a} <=> $href->{$i}{$b}} keys %{$href->{$i}})[0];
		
		my @bases = keys %{$href->{$i}};
		my $commonBase = '';
		if ( scalar @bases > 1 ) {
			$commonBase = join("", ("[", %{$href->{$i}}, "]"));
		} else {
			$commonBase = $bases[0]; # Should be only one base!
		}
		
		push(@cons, $commonBase);		
		
	}
	
	my $c = join("", @cons);
	
	return $c;

} # end sub computeConsensusLong


sub identifyOffTarget {
	my $seq = shift @_;
	my $logger = Log::Log4perl->get_logger('debarcer');

	my ($bc_position, $barcode) = &extractBarcode($seq);
	my $whichAmplicon = 'UNKNOWN';
	my $ampliconSeq = '';
	
	my $primerList = '';
	
	foreach my $ampID ( keys %primerSets ) {
		my ($fwd, $rev) = @{$primerSets{$ampID}};
		
		if ( $seq =~ m/$fwd/g ) {
			$primerList .= "$ampID" . '_F|' . pos($seq) . ';';
			pos($seq) = 0;
		}
		if ( $seq =~ m/$rev/g ) {
			$primerList .= "$ampID" . '_R|' . pos($seq) . ';';
			pos($seq) = 0;
		}
	
		# Now do the antisense FIXME: what were we doing before??
		$fwd = &revcom($fwd);
		if ( $seq =~ m/$fwd/g ) {
			$primerList .= "as$ampID" . '_F|' . pos($seq) . ';';
			pos($seq) = 0;
		}
		$rev = &revcom($rev);
		if ( $seq =~ m/$rev/g ) {
			$primerList .= "as$ampID" . '_R|' . pos($seq) . ';';
			pos($seq) = 0;
		}
	}
	
	$logger->debug("\n$seq\nAMPLICON= $whichAmplicon $barcode $primerList\n");
	
	return ($whichAmplicon, $bc_position, $barcode, $primerList);	
}



sub listPrimers {
	my $aref = shift @_;  # This is @ampliconsFile - FIXME what is, the input parameter(s)?
	
	my %p = ();
	foreach ( @$aref ) {
		next if (/^#/); chomp;
		my ( $n, $id, $f, $r, $seq, $note ) = split /\t/;
		$r = &revcom($r);
		@{$p{$id}} = ( $f, $r );
	}
			
	return %p;
}

sub ampliconLengths {
	my $aref = shift @_;  # This is @ampliconsFile - FIXME same question as above
	
	my %h = ();
	foreach ( @$aref ) {
		next if (/^#/); chomp;
		my ( $n, $id, $f, $r, $seq, $note ) = split /\t/;
		$h{$id} = length($seq);
	}
	return %h;
	
} # end sub listAmplicons


sub ampliconSequences {
	my $aref = shift @_;  # This is @ampliconsFile - FIXME same as above
	
	my %h = ();
	foreach ( @$aref ) {
		next if (/^#/); chomp;
		my ( $n, $id, $f, $r, $seq, $note ) = split /\t/;
		$h{$id} = $seq;
	}
	return %h;
	
} # end sub ampliconSequences


sub extractBarcodeQuick {
	
=pod

extractBarcodeQuick is a fast way to extract barcodes from bam file alignments
and assumes that the UID is present at either the start or end of the sequence
depending on whether the alignment is on the forward or reverse strand.

If the UID is not present at position 0, no attempt to salvage the read is made.

=cut
	my $logger = Log::Log4perl->get_logger('debarcer');

	my $seq = shift @_;
	my $revseq = &revcom($seq);
	my $bc_pos = '';
	my $barcode = '';
	my $isRev = 0;
	
	my $barcodeSuffix = "ATGGGAAAGAGTGT"; # Modification to barcodeSuffix to facilitate fuzzy matching.  Do not use regex with length statement below.
	
	# my $seqHead = substr($seq, 0, 26); # This is 12 + 8 + 6 = UID + spacer + GAGTGT motif
	# but this is faster
	if ( distance(substr($seq, 12, 14), $barcodeSuffix) <= 2 ) { # Check for the motif
		$barcode = substr($seq, 0, 12);
		$bc_pos = 0;
	} elsif ( distance(substr($revseq, 12, 14), $barcodeSuffix) <= 2 ) {
		$barcode = substr($revseq, 0, 12);
		$bc_pos = 0;
		$isRev = 1;
	}
	
	$logger->debug("$bc_pos $barcode $isRev $seq\n");
	return ($bc_pos, $barcode, $isRev);
	
} # End sub extractBarcodeQuick

sub extractBarcode {

=pod

extractBarcode was the original way of finding barcodes with SaferSeq amplicons.
Deprecated.  Use extractBarcodeQuick instead.

=cut

	my $seq = shift @_;
	my $logger = Log::Log4perl->get_logger('debarcer');

	# FIXME whoa - what happened here? should the barcodeSuffix be in the config file or command-line parameters?

	# my $barcodeSuffix = "ATGGGAAAGAGTGTCC"; # From Hair-Seq deck June 18th.
	# my $barcodeSuffix = "ATGGGAAAGAGTGTGG"; # From Hair-Seq deck June 18th.
	# my $barcodeSuffix = "ATGGGAAAGAGTGT[G{2}|C{2}]"; # From Hair-Seq deck June 18th.
	# my $barcodeSuffix = "ATGGGAAAGAGTGT[GC]{2}"; # From Hair-Seq deck June 18th.  Correction to regex
	my $barcodeSuffix = "ATGGGAAAGAGTGT"; # Modification to barcodeSuffix to facilitate fuzzy matching.  Do not use regex with length statement below.
	# my $barcodeSuffix = "ATGGGAAAGA"; # Modification to barcodeSuffix to facilitate fuzzy matching.  Do not use regex with length statement below.
	
	my $bcSuffixMatch = &fuzzyMatch2($barcodeSuffix, $seq);
	# my $bcSuffixMatch = $barcodeSuffix;
	
	$seq =~ m/($bcSuffixMatch)/g;
	my $bc_pos = (pos $seq) - 12 - length($1); # pos - barcodeLength - barcodeSuffixLength
	# FIXME: is this a logging call or an error report?
	# print STDERR "BarcodeSuffixFound: $1 $bc_pos\n";
	my $barcode = substr( $seq, $bc_pos, 12);
	
	# FIXME: is this a logging call or an error report?
	# print "\nSEQ= $seq\n" . "BCS= $barcodeSuffix\n" . "BCP= $bc_pos\n" . "BC=  $barcode\n";
	
	return ($bc_pos, $barcode);
}


sub detectAmplicons {

	my ($bamfile, $plexity) = @_;
	my %allAmplicons = ();
	my $logger = Log::Log4perl->get_logger('debarcer');
	
	open BAM, "samtools view -s 0.01 -X $bamfile | cut -f 3,4 |"; # Sample a fraction of the file.
	while (<BAM>) {
		
		chomp;
		my @line = split /\t/;
		my $chrom = $line[0];
		next if ( $chrom eq "*" );
		my $chromStart = $line[1];
		
		my $ampName = "$chrom:$chromStart";
		$allAmplicons{$ampName}++;
		
	}
	close BAM;
	
	my @validKeys = sort { $allAmplicons{$b} <=> $allAmplicons{$a} } keys %allAmplicons;
	my %validAmplicons = ();
	$logger->info("Setting valid amplicons based on plexity = $plexity:\n");
	for (my $i = 0; $i < $plexity; $i++) {
		$logger->info("$validKeys[$i]\n");
		$validAmplicons{$validKeys[$i]}++;
	}
		
	return %validAmplicons;
		
}


sub getPositionAliases {

=pod

=head2 getPositionAliases

Read in a table of position-amplicon names useful for making output more human-readable

=cut

	my $infile = shift @_;  # Usually ./amplicon_tables/all_amplicons.txt
	my %h = ();
	my $logger = Log::Log4perl->get_logger('debarcer');
	
	my @header = '';
	
	open INFILE, $infile;
	while (<INFILE> ) {
		chomp;
		
		if ( /^#/ ) { # collect header
			$_ =~ s/^#//;
			@header = split /\t/;
			next;
		}

		my %line = ();
		@line{@header} = split /\t/;

		if (  $line{"AmpliconStart"} ne "" ) {
			$h{$line{"AmpliconStart"}} = $line{"AmpliconName"};
		}

	}
	close INFILE;
	
	return %h;
}

sub loadAmpliconData {

=pod

=head2 loadAmpliconData

Read in a table of information pertaining to the amplicons being analyzed

=cut

	my $infile = shift @_;  # Usually ./amplicon_tables/all_amplicons.txt
	my $href = shift @_;

	my @header = '';

	open INFILE, $infile;
	while (<INFILE> ) {
		chomp;
		
		if ( /^#/ ) { # collect header
			$_ =~ s/^#//;
			@header = split /\t/;
			next;
		}

		my %line = ();
		@line{@header} = split /\t/;

		if (  $line{"TargetWindow"} ne "" ) {
			my ($chrom, $start, $end) = split(/[:-]/,  $line{"TargetWindow"});
			$href->{$line{"AmpliconName"}}{"TargetWindow"}{"chrom"} = $chrom;
			$href->{$line{"AmpliconName"}}{"TargetWindow"}{"start"} = $start;
			$href->{$line{"AmpliconName"}}{"TargetWindow"}{"end"} = $end;
		}

	}
	close INFILE;
}

=pod 

=head2 trimSeq
 
Remove the end of the sequence with poor QUAL values
 
=cut

sub trimSeq {
	my ($seq, $quals) = @_;
	$quals =~ m/(#+$)/g;
	my $badQualString = $1;
	$seq = substr($seq, 0, length($seq) - length($badQualString) + 1);
	return $seq;
}

=head2 unique
 
Return unique items from an array
 
=cut

sub unique {
	my %t = ();
	foreach (@_) { $t{$_}++; }
	return keys %t;
}

sub revcom {
	my $seq = shift;
	$seq =~ tr/ACTGN/TGACN/;
	$seq = reverse $seq;
	return $seq;
} # end sub revcom

1;
