#!/usr/bin/perl
use strict;
use Data::Dumper;
(open my $TTY,"/dev/tty") || die "unable to open keyboard input"; ## will allow code to be modified with <$TTY>, to hold for keyboard input.  <STDIN> will not work s this script reads from a stream

=pod

Identify barcodes that are probably mutated UIDs

=cut

my $levenshtien_binary = "$ENV{'BHOME'}/bin/levenshtien_stream";

# Test the Levenshtein function
#print distance("foo","four");
# prints "2"

my @UIDs = '';

my %prefix_hash = ();
my %suffix_hash = ();
my %depth_hash = ();

my @inputstream = <>;   # Ugly
my @amplicons = &getAmplicons(\@inputstream);
print STDERR "Finding edits in @amplicons\n";

## LEH : collects ALL 6mer prefix and suffix for all the uids associated with an amplicon
## LEH : does this account for the correct barcode length for suffixes?  concern that this is tuned to SimSenSeq and not HaloPlex
# my $amplicon = "SMA548";
foreach ( @inputstream ) {
	chomp;
	my ($thisAmplicon, $UID, $depth, $seq, $seqMarkup) = split /\t/;
	next if ( $thisAmplicon eq "UNKNOWN" );
	
	# 	push(@UIDs, $UID);
	my $UID_pre = substr($UID, 0, 6);
	### my $UID_suff = substr($UID, 7, );  this is wrong, i think the intention is to capture the last 6 bases, probably desinged aroudn 12mers.  haloplex is 8 mers
	my $UID_suff=substr($UID,length($UID)-6,6);
	push(@{$prefix_hash{$thisAmplicon}->{$UID_pre}}, $UID);
	push(@{$suffix_hash{$thisAmplicon}->{$UID_suff}}, $UID);
	$depth_hash{$thisAmplicon}->{$UID} = $depth;
}


foreach my $thisAmplicon ( keys %depth_hash ) {
	
	print STDERR "Finding edits in $thisAmplicon\n";
        ### LEH :goes through all prefixes and sufficxes associate with an amplicon identifies those that are close (1bp?)
        ### LEH : uses the levenstien binary
	foreach my $href ( ( \%prefix_hash, \%suffix_hash ) ) {  ### goes through both hashes
		foreach my $prefix ( keys %{$href->{$thisAmplicon}} ) {
			my @UIDs = @{$href->{$thisAmplicon}{$prefix}};
			
			my $nUIDs = scalar @UIDs;
			next if ( $nUIDs == 1 );
			# print "$prefix: $nUIDs\n";
			my $UIDs_out = join("\n", @UIDs);
			my @closeUIDs = `echo "$UIDs_out" | $levenshtien_binary`;
			
			
			
			next unless ( scalar(@closeUIDs) );

			print "$prefix\n";
			print Dumper(@UIDs);<$TTY>;
			print Dumper(@closeUIDs);<$TTY>;
			foreach ( @closeUIDs ) {
				next unless (/^Edit/);
				chomp;
				my ( $descriptor, $bc1, $bc2 ) = split("\t", $_);
				my $d1 = $depth_hash{$thisAmplicon}{$bc1};
				my $d2 = $depth_hash{$thisAmplicon}{$bc2};
			        ### LEH : for each pair, compares depth
				### LEH : masks	(not corrects) barcode if ratio of depths is <0.2 and depth >5
				my $ratio;
				if ( $d2 >= $d1 ) {
					$ratio = $d1 / $d2;
					my $mask = 'OK';
					if ( ($ratio < 0.2) & ($d2 > 5) ) { # This is the masking heuristic
						$mask = "Mask";
						} elsif ( $d1 == 1 ) {
						$mask = "Mask";
						}
					printf("DuplicateIn:\t%s\t%s\t%d\t%s\t%d\t%.3f\t%s\n", 
						$thisAmplicon, $bc1, $d1, $bc2, $d2, $ratio, $mask);
					} else {
					$ratio = $d2 / $d1;
					my $mask = 'OK';
					if ( ($ratio < 0.2) & ($d1 > 5) ) { # This is the masking heuristic
						$mask = "Mask";
						} elsif ( $d2 == 1 ) {
						$mask = "Mask";
						}
					printf("DuplicateIn:\t%s\t%s\t%d\t%s\t%d\t%.3f\t%s\n", 
						$thisAmplicon, $bc2, $d2, $bc1, $d1, $ratio, $mask);
				}
			}
		}
	}
}
		
exit;

###########################################

sub getAmplicons {
	my $aref = shift @_;
	my %ampHash = ();
	foreach ( @$aref ) {
		chomp;
		my ($amp, $uid, $depth, $seq, $seqMarkup) = split /\t/;
		$ampHash{$amp}++;
		}
	my @amps = keys %ampHash;
	return @amps;
}

