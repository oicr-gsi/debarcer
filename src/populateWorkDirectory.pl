#!/usr/bin/perl
use strict;

=pod

This script generates the 'launchDebarcer.sh' and 'relaunchAllDebarcers.sh' scripts used to 
kick off a Debarcer run

Usage:

  perl $BHOME/src/populateWorkDirectory.pl <directory with the fastqs to analyze>


(c) Paul Krzyzanowski 2015-2016
Email: paul.krzyzanowski@oicr.on.ca

=cut

use Getopt::Long;
use File::Basename;
use Cwd;

my $fastq_source = shift @ARGV or die "\nNeed to provide a directory with fastq.gz files to process...\n\n";
my @fastqs = `ls $fastq_source/*_R1_001.fastq.gz`;

# Set plexity if given
my $plexity = shift @ARGV;
$plexity = ( $plexity ) ? $plexity : 50;
print STDERR "Using plexity $plexity\n";

my $module_version = "dev";
&GetOptions(
	"module=s" => $module_version
);


open RELAUNCHSCRIPT, ">relaunchAllDebarcers.sh";
print RELAUNCHSCRIPT "#!/bin/bash\n\n";
print RELAUNCHSCRIPT "LAUNCHDIR=`pwd`\n\n";


foreach my $fastq ( @fastqs ) {

	chomp $fastq;
	my $thisdir = cwd();
	my $fastq_abs_name = $thisdir . "/" . $fastq;

	my $samplename = basename($fastq, "_R1_001.fastq.gz");
	my $basefile = basename($fastq);
	print STDERR "$samplename $basefile\n";

	system("mkdir -p results/$samplename");
	
	open LAUNCHSCRIPT, ">results/$samplename/launchDebarcer.sh";
	print LAUNCHSCRIPT "qsub -N \"Dbarc_$samplename\" -l h_vmem=32G -cwd -b y \"module load debarcer/$module_version; runDebarcer.sh -r -f $fastq_abs_name -n $samplename -o .\"";
	close LAUNCHSCRIPT;
	
	open CONFIGSCRIPT, ">results/$samplename/debarcer.conf";

	print CONFIGSCRIPT "# Debarcer override file.  The master process stores these variables in the config array.\n";
	print CONFIGSCRIPT "plexity=$plexity\n";

	close CONFIGSCRIPT;
	
	print RELAUNCHSCRIPT <<BLOCK
cd results/$samplename/
echo results/$samplename/
rm Dbarc*; ./launchDebarcer.sh
cd -

BLOCK
	;

}
	
close RELAUNCHSCRIPT;

exit;
