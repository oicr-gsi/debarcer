#!/usr/bin/perl
use strict;

=pod

Combine HaloplexHS barcode fastqs with read fastqs

  perl combineHaloplexFastqs.pl --R1file [] --R2file [] --barcodefastq [] --outputdir []
  
  --R1file is the fastq.gz that contains your actual Read 1s
  --R2file is the fastq.gz that contains your actual Read 2s
  --barcodefastq is the fastq.gz containing 'reads' that are the HaloplexHS UMI

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
	"outputdir=s"		=> \$args{"outputdir"},
	"R1file=s"			=> \$args{"R1file"},
	"R2file=s"			=> \$args{"R2file"},   
	"barcodefastq=s"	=> \$args{"barcodeFq"}   
);




open my $BARCODES, "gunzip -c ".$args{"barcodeFq"}." |";
open my $R1IN, "gunzip -c ".$args{"R1file"}." |";
open my $R2IN, "gunzip -c ".$args{"R2file"}." |";


my %FHOUT;
my $R1outputfile = $args{"outputdir"} ."/". basename($args{"R1file"});
$R1outputfile =~ s/_R\d_/_R1_/;  # Ensure the file indicated R1
(open $FHOUT{R1}, "| gzip -c > ". $R1outputfile) || die "unable to open $R1outputfile outputfile";


my $R2outputfile = $args{"outputdir"}."/". basename($args{"R2file"});
$R2outputfile =~ s/_R\d_/_R2_/;  # Ensure the file indicated R2
(open $FHOUT{R2}, "| gzip -c > ". $R2outputfile) || die "unable to open $R1outputfile outputfile";


my $reccount=0;
until(eof($BARCODES)){
	$reccount++;
	my %rec;
	%{$rec{bc}}=get_fastq_rec($BARCODES);
	%{$rec{R1}}=get_fastq_rec($R1IN);
	%{$rec{R2}}=get_fastq_rec($R2IN);
	

	### TO DO : Add code to check concordance of records
	###       should be a matter of the first part of the header matching in all fastq files
	
	my $barcode=$rec{bc}{read};
	for my $R(qw/R1 R2/){
		my @header=split /\s/,$rec{$R}{head};
		$header[0].=":HaloplexHS-".$barcode;
		$rec{$R}{head}=join(" ",@header);
		print {$FHOUT{$R}} join("\n",@{$rec{$R}}{qw/head read head2 qual/}) . "\n";
	}
}
print STDERR "$reccount fastq records processed\n";



sub get_fastq_rec{
	my ($FH)=@_;
	my %rec;
	for my $line(qw/head read head2 qual/){
		$rec{$line}=<$FH>;
		chomp($rec{$line});
	}
	return %rec;
}







	
	
	
	


 
