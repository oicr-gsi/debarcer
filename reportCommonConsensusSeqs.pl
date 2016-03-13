

=pod

Summarize the stream from a command like

  perl classifyReadsToAmplicons.pl -fastqgz=../fastqs/Sample_01.R1.fastq.gz | perl reportCommonConsensusSeqs.pl
  
  perl classifyReadsToAmplicons.pl -fastqgz=../fastqs/Sample_01.R1.fastq.gz | perl reportCommonConsensusSeqs.pl --familydepth=20
  
  
  
=cut


use Getopt::Long;

my %args = ();
GetOptions(
	"familydepth=s" => \$args{"familydepth"}
);

my $family_depth = ( $args{"familydepth"} ) ? $args{"familydepth"} : 3;

my %Amps = ();

while (<>) {

	chomp;
	my ($amplicon, $barcode, $count, $consensus, $consensusLong) = split /\t/;
	
	next if ( $count < $family_depth );
	
	$Amps{$amplicon}->{$consensus}++;

}

foreach my $amplicon ( keys %Amps ) {
	my $counter = 15;
	foreach my $seq ( sort {$Amps{$amplicon}{$b} <=> $Amps{$amplicon}{$a}} keys %{$Amps{$amplicon}} ) {
		printf("%s\t%s\t%s\n", $amplicon, $seq, $Amps{$amplicon}->{$seq} );
		$counter--;
		#last if ( $counter == 0 );
	}
	
}

exit;
	
