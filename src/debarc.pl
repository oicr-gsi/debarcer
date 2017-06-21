##!/usr/bin/perl
use strict;

=pod

=head1 generateConsensusFromBAM.pl

Description

=head2 Usage

   perl generateConsensusFromBAM.pl --bam=
   	   
=head2 Author

Paul M Krzyzanowski
pmkrzyzanowski@gmail.com
(c) 2014-2016

=cut

# Setup general modules
# use Cwd 'abs_path';
# use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use Debarcer;
use Getopt::Long;
use Data::Dumper;
use Bio::DB::Sam;
use JSON::XS qw(encode_json decode_json);
use Config::General qw(ParseConfig);

my $DBROOT = "$FindBin::Bin/../"; # Because script is in $DBROOT/src/


my %opts = ();

GetOptions(
	"config=s" 		=> \$opts{"configfile"},
	"bam=s" 		=> \$opts{"bam"},
	"sampleID=s" 	=> \$opts{"sampleID"},
	"consDepth=s" 	=> \$opts{"consDepth"},  
	"plexity=s" 	=> \$opts{"plexity"},
	"offset=i"		=> \$opts{"offset"},  
	# "strictCons" 	=> \$opts{"strictCons"},  # Strict consensus
	"downsample=s" 	=> \$opts{"downsample"},
	"justUIDdepths"	=> \$opts{"justUIDdepths"},
	"justTargets" 	=> \$opts{"justTargets"},
	"test"			=> \$opts{"test"},  # test mode
	"sites=s"		=> \$opts{"sitesfile"},
	"site=s"		=> \$opts{"site"},   ## a signle site to process
	"UIDdepths"		=> \$opts{"UIDdepths"},
	"basecalls" 	=> \$opts{"basecalls"},
	"out=s"			=> \$opts{output_folder},
	"familysizes=s" => \$opts{familysizes},
	"help" 			=> \$opts{help},
);

### configuration options are overwritten by command line options
my %config = ParseConfig($opts{"configfile"});
map{$opts{$_} = $opts{$_} || $config{$_}} keys %config;
validate_options(\%opts);


print STDERR "--- Starting debarc.pl ---\n";

# Section to load parameters from a Config::Simple format file
my $nSites = $opts{"plexity"};
### set up subfolders in the output folder
my $table_folder=$opts{output_folder} ."/tables";
mkdir($table_folder) unless(-d $table_folder);

### sets the amplicon table for loading a list of expected amplicons for the library
### sets it to a default value if not provided in the configuration
###  BETTER to not use this variable, but instead change the config parameter, or set this as an argument. this is a useless variable
my $ampliconTable = ( $config{"ampliconTable"} ) ? $config{"ampliconTable"} : "$DBROOT/data/all_amplicons.txt" ;

$opts{"justTargets"} = 1;   ### this is not currently functional, code is commented out.  need to revies
my $consensusDepth = $opts{"consDepth"};
print STDERR "Using Consensus Depth = $consensusDepth and plexity = $nSites\n";

#### CONSTANTS
my @SNVtypes=qw/A C G T D I N/;  ### IS THIS ALL TYPES?

### globals
my $inputSeqCount = 0;
my $familySitesSeqCount = 0;


my $prefix="ALL";
my %sites;
if($opts{site}){
	$sites{$opts{site}}=1;
	($prefix=$opts{site})=~s/:/_/;
	
}else{
	### the list of sites should be provided in a file.
	### if the file does not exist, then it needs to be created
	%sites=load_sites($opts{"sitesfile"},$nSites,$opts{"bam"});
	
}
my $count=scalar keys %sites;
print STDERR $count ." family sites loaded\n";




## if these flags aren't set, then proceed no further
exit unless ($opts{UIDdepths} || $opts{basecalls});

### alisases for annotation, this will identify the position of each of the amplicon positions in the genome and the related gene
print STDERR "loading site aliases from $ampliconTable\n";
my %siteAliasTable = &Debarcer::getPositionAliases($ampliconTable);

#### This is a hash of genes and "target windows.  not sure what this is used for???/
### Appears to be used for the justTarget mode, which is currently disabled.  review later
print STDERR "loading amplicon info from $ampliconTable\n";
my %ampliconInfo = ();
&Debarcer::loadAmpliconData($ampliconTable, \%ampliconInfo);

#### I DON"T UNDERSTAND THIS"
print STDERR "identifying invalid barcodes\n";
my %invalidBarcodes = &loadBarcodeMaskFile($opts{"sampleID"});
#print Dumper(\%invalidBarcodes);<STDIN>;



### this does not appear to be used?
#my $bam = $sam->bam();
#my $header = $bam->header();

### open files as needed
my ($UIDDEPTHFILE,$CONSENSUSFILE,$POSITIONFILE);
if ($opts{"UIDdepths"}){
	my $uidDepthFile = $table_folder. "/" . $prefix . "." . $opts{"sampleID"} . ".UIDdepths.txt.gz";
	open $UIDDEPTHFILE, "| gzip -c > $uidDepthFile";
	print STDERR "opening $uidDepthFile\n";
}

if ($opts{"basecalls"}){
	my $ConsensusFile = $table_folder. "/" . $prefix . "." . $opts{"sampleID"} . ".consensusSequences.txt.gz";
	open $CONSENSUSFILE, "| gzip -c > $ConsensusFile";
	
	my $positionFile = $table_folder. "/" . $prefix . "." . $opts{"sampleID"} . ".positions.txt";
	open $POSITIONFILE, ">",$positionFile;
	print $POSITIONFILE join("\t", "#Chrom", "Alias", "Position", "Ref","AmpPosition","FSize",join("\t", @SNVtypes), "Depth","RefProp")."\n";
	

}

print STDERR "begin parsing bam file $opts{bam}\n";
my $sam = Bio::DB::Sam->new(-bam => $opts{bam}, -fasta => $config{refgenome} );

my $sitecount=0;
my @AmpliconIDs=sort {$sites{$b}<=>$sites{$a}} keys %sites;
@AmpliconIDs=@AmpliconIDs[0..($nSites-1)] if($nSites>scalar @AmpliconIDs);

for my $AmpliconID(@AmpliconIDs){
	$sitecount++;
	my($chrom,$start,$stop,$count)=split /[:-]/,$AmpliconID;
	print "extracting reads at $AmpliconID\n";

	### sitedata will store the results of parsing the bam file for reads at the specific site
	### it will include uid family stats, consensus information and base calls by positions
	#my %sitedata=get_site_data($chrom,$start,$stop,$sam,$AmpliconID,$opts{offset});
	my %sitedata=get_site_data_pe($chrom,$start,$stop,$sam,$AmpliconID,$opts{offset});

	
	
	
	
	#collapse_uids($sitedata{uids});
	
	my $ampliconReportingThreshold = &calculateAmpliconReportingThreshold_rev(%{$sitedata{basecounts}});
	my $uidcount=scalar keys %{$sitedata{uids}};
	print "$uidcount UIDs found at this site\n";
	
	( my $printableAmpID = $AmpliconID ) =~ s/\t/:/;
	
	
	
	
	# print to the UID depth file if the UIDdepths flag is set.
	print_uid_depths($UIDDEPTHFILE,$printableAmpID,$sitedata{uids}) if($opts{UIDdepths});
	

	next unless ($opts{basecalls});  ### proceed no further if basecalls are not asked for
	
	### pass the sitedata as a reference.  it will be modified by this function to store consensus sequences
	
	generate_consensus_data($AmpliconID,$CONSENSUSFILE,\%sitedata,$opts{familysizes});
	
	
	#print Dumper(%sitedata);exit;
	#my $ampliconName = $siteAliasTable{$AmpliconID};
	#unless ($ampliconName) {
	#	$ampliconName = &Debarcer::generateAmpliconName($AmpliconID);
	#}
	### BETTER WRITTEN AS
	my $ampliconName = $siteAliasTable{$AmpliconID} || &Debarcer::generateAmpliconName($AmpliconID);

	my $probableRefBase="N";  ## this shoudl be the reference base
	for my $pos( sort{$a<=>$b} keys %{$sitedata{basecalls}}){

		my $refpos=$start+$pos;
		my $refbase=$sam->seq($chrom,$refpos,$refpos);
		
		

		my @counts=map{ $sitedata{basecalls}{$pos}{raw}{$_} || 0 } @SNVtypes;
		my $rawdepth=$sitedata{basecalls}{$pos}{raw}{depth};
		
		
		my $refcount=$sitedata{basecalls}{$pos}{raw}{$refbase} || 0;
		my $reffreq=$rawdepth ? $refcount/$rawdepth : 0;

		printf $POSITIONFILE ("%s\t%s\t%d\t%s\t%d\t%s\t", $chrom, $ampliconName, $refpos, $refbase,$pos,"raw");
		print $POSITIONFILE join("\t", @counts) . "\t$rawdepth\t$reffreq\n";

		#print "$pos\n";
		#print Dumper($sitedata{basecalls}{$pos});<STDIN>;


		
		my @uidlevels=split /,/,$opts{familysizes};
		
		for my $uidlevel(@uidlevels){
			
			
				my @counts=map{ $sitedata{basecalls}{$pos}{consensus}{$uidlevel}{$_} || 0 } @SNVtypes;
				my $consdepth=$sitedata{basecalls}{$pos}{consensus}{$uidlevel}{depth};
				
				my $refcount=$sitedata{basecalls}{$pos}{consensus}{$uidlevel}{$refbase} || 0;
				my $reffreq=$consdepth ? $refcount/$consdepth : 0;				

				printf $POSITIONFILE ("%s\t%s\t%d\t%s\t%d\t%s\t", $chrom, $ampliconName, $refpos, $refbase,$pos,"cons_$uidlevel");
				print $POSITIONFILE join("\t", @counts) . "\t$consdepth\t$reffreq\n"

		}

	}
	
	
	next;
	
	#######
	
	
	
	
	my %calltable=get_position_calls($AmpliconID,\%sitedata,$ampliconReportingThreshold);
		
	
	for my $pos ( sort {$a <=> $b} keys %calltable){
		
		my ($rawDepth, $consDepth) = ( 0, 0);
		my ($probableRefBase, $n) = ('', 0);
		for my $base ( @SNVtypes ) {
			my $snvRawDepth = $calltable{$pos}{"raw"}{$base};
			$rawDepth += $snvRawDepth; # Increment rawDepth
			if ( $snvRawDepth > $n ) { # If this is the highest raw depth observed, make the base the probableBase
				$probableRefBase = $base;
				$n = $snvRawDepth;
			}
			
			$consDepth += $calltable{$pos}{"consensus"}{$base}; # Increment consDepth
		}
		
		printf $POSITIONFILE ("%s\t%s\t%s\t%s", $chrom, $ampliconName, $pos, $probableRefBase);
		printf $POSITIONFILE ("\t%d", $calltable{$pos}{"raw"}{$_}) foreach ( @SNVtypes );
		printf $POSITIONFILE ("\t%d", $rawDepth);
		printf $POSITIONFILE ("\t%d", $calltable{$pos}{"consensus"}{$_}) foreach ( @SNVtypes );
		printf $POSITIONFILE ("\t%d", $consDepth);
		print $POSITIONFILE "\n";
	}
	#last if $sitecount>2;
}

print STDERR "Raw reads in family sites read from $opts{bam}: $familySitesSeqCount\n";



print STDERR "--- generateConsensusFromBAM.pl Complete ---\n";

close $UIDDEPTHFILE if ($opts{"UIDdepths"});
close $CONSENSUSFILE if ($opts{"basecalls"});
close $POSITIONFILE  if ($opts{"basecalls"});

exit;

#################################################################################################



sub load_sites{
	### called to either load a list of sites that has already been generated, 
	### or a call to the function that will generate this list and save to the file
	
	my ($file,$nSites,$bamfile)=@_;
	my %sites;
	### if the file exists, then read in the sites
	if(-e $file){
		(open my $FSITES,"<",$file) || die "unable to open sites file";
		while(<$FSITES>){
			chomp;
			my($site,$size)=split /\t/;
			$sites{$site}=$size;
		}
		close $FSITES;
	}else{
		print STDERR "sites file not found. Identifying Family sites : $nSites\n";
		#(open my $FSITES,">",$opts{"sampleID"}.".familysites") || die "unable to open family sites";
		%sites = &identifyFamilySites_PE($bamfile, $nSites);

		#### write the sties to a file
		(open my $FSITES,">",$file) || die "unable to open sites file";
		## sort by size
		for my $site(sort { $sites{$b} <=> $sites{$a} } keys %sites){
			print $FSITES "$site\t$sites{$site}\n";
		}
		close $FSITES;
	}
	return %sites;
}

sub identifyFamilySites {
	
	### will sample from a bam file to identify the most frequenly seen sites, printing out the top N to a file
	my $inBam = shift @_;
	my $nSites = shift @_;
	my %familySites = ();

	print STDERR $ENV{"SAMTOOLSROOT"}."\n";
	my $SAMTOOLSBINARY = $ENV{"SAMTOOLSROOT"}."/bin/samtools";
	
	
	#### better approach since already using bash tools
	#my @allSites = `$SAMTOOLSBINARY view -s 0.1 $inBam | cut -f 3,4 | sort | uniq -c | sort -nr`;
	### this should be faster.  counts already done and sorted
	
	my @allSites = `$SAMTOOLSBINARY view -s 0.1 $inBam | cut -f 3,4`;
	foreach my $site ( @allSites ) { 
		chomp $site; 
		$site =~ s/\t/:/; $familySites{$site}++ unless($site=~/\*/); 
	}
	my @goodSites = (sort { $familySites{$b} <=> $familySites{$a} } keys %familySites);
	@goodSites = @goodSites[0..$nSites];  # Take the top n-1
	print STDERR "Compiling info for Family Sites:\n";
	print STDERR "Note:  If present, the '* 0' (i.e. unmapped) site has been dropped\n";
	
	my %goodHash = ();
	@goodHash{@goodSites} = @familySites{@goodSites};
	return %goodHash;
	
}

sub identifyFamilySites_PE {
	
	### will sample from a bam file to identify the most frequenly seen sites, printing out the top N to a file
	my $inBam = shift @_;
	my $nSites = shift @_;
	my %familySites = ();

	print STDERR $ENV{"SAMTOOLSROOT"}."\n";
	my $SAMTOOLSBINARY = $ENV{"SAMTOOLSROOT"}."/bin/samtools";
	
	
	#### better approach since already using bash tools
	#my @allSites = `$SAMTOOLSBINARY view -s 0.1 $inBam | cut -f 3,4 | sort | uniq -c | sort -nr`;
	### this should be faster.  counts already done and sorted
	
	#my @allSites = `$SAMTOOLSBINARY view -s 0.1 $inBam | cut -f 3,4,7,9`;   ### chrom
	
	## DON"T SAMPLE< AND ONLY GET READ1 f2=proper pair, F16=NOT on the reverse strandr
	(open my $BAM,"$SAMTOOLSBINARY view -f 2 -F 16 $inBam | cut -f 2,3,4,7,9 |") || die "could not open bam stream";
	while(my $site=<$BAM>){
		chomp $site; 
		my($flag,$rname,$pos,$rnext,$tlen)=split /\t/,$site;
		next unless($rnext eq "=");
		my $pos2=$pos+$tlen-1;
		my $siteid="${rname}:${pos}-${pos2}";
		$site =~ s/\t/:/; $familySites{$siteid}++;
	}
	close $BAM;
	my @goodSites = (sort { $familySites{$b} <=> $familySites{$a} } keys %familySites);
	@goodSites = @goodSites[0..$nSites];  # Take the top n-1
	print STDERR "Compiling info for Family Sites:\n";
	#print STDERR "Note:  If present, the '* 0' (i.e. unmapped) site has been dropped\n";
	
	my %goodHash = ();
	@goodHash{@goodSites} = @familySites{@goodSites};
	return %goodHash;
	
}


sub get_site_data{
	
	### extracts read at a specific start point and stores data about these reads inot thhe has
	
	my($chrom,$start,$sam,$AmpliconID,$offset)=@_;
	
	my %data;
	## GET ALL ALIGNMENTS THAT OVERLAP THIS that cover this positions
	my $segment = $sam->segment($chrom,$start,$start);
	#my @alignments = $segment->features;
	## LIMIT TO ALIGNMENTS AT THIS EXACT START POSITION	
    my @alignments = $segment->features(-filter => sub { my $alignment = shift;return $alignment->start==$start;});
	my $count_exact=scalar @alignments;
	my $count_offset=0;
	if($offset){
		my @offset_alignments=$segment->features(-filter=>sub{ 
			my $alignment=shift;
			my $lower=$start-$offset;
			my $upper=$start+$offset;
			return $alignment if( ($alignment->start>=$lower) && ($alignment->start<=$upper) && ($alignment->start!=$start) );
		});
		$count_offset=scalar @offset_alignments;
		push(@alignments,@offset_alignments);
	}

	### DEBUG, REMOVE for production
	#my $max=10000;
	#my $count=scalar @alignments;
	#$max=$count<$max ? $count : $max;
	#@alignments=@alignments[0...$max];	
	
	
	my $count_all=scalar @alignments;
	print STDERR "$count_all alignments for this site, $count_exact + $count_offset within $offset bases\n";
	
	### PROCESS EACH AlIGNMENT
	#my %align_lengths;
	for my $alignment(@alignments){
			
		#print "$inputSeqCount $chrom $chromStart";<STDIN>;
		$data{SitesSeqCount}++;
		

		
		my $barcode = '';
		my $bc_position = 0;
		my $read_name = $alignment->query->name();
		if ($read_name =~ /HaloplexHS/) { # We have a HaloplexHS read
			$read_name =~ m/HaloplexHS-([ACTG]{10})/;
			$barcode = $1;
		} else { # We have a SiMSenSeq read
			my $read_dna = $alignment->query->dna();
			($bc_position, $barcode) = &Debarcer::extractBarcodeQuick($read_dna);
		}
		next unless ( $bc_position == 0 );
		# Skip the barcode if it's in the mask hash loaded from the maskfile
		if ( exists $invalidBarcodes{$AmpliconID}->{$barcode} ) {
			# print "Skipping invalidBarcodes{$AmpliconID}->{$barcode}\n"; 
			next unless ( $opts{"justUIDdepths"} ); # Do not skip an invalid barcode if we're only generating the UID depth file.
		}
	
		# Count this instance of the barcode family
		$data{uids}{$barcode}{count}++;
	    
		### if not askig for basecalls, tehn no point going any futher
		next unless($opts{basecalls});
		
		my $ref_dna   = $alignment->dna; 
		my $query_dna = $alignment->query->dna; # query sequence bases
		#print "$query_dna\n";
		#print "$ref_dna\n";<STDIN>; 
		#$align_lengths{length($ref_dna)}++;
	
		# Determine the base call or insertion/deletion status wrt to the start of the alignment, i.e.
		# chromosomal position given in $chromStart.  Therefore $basecalls[0] is for $chromStart + 0
		my @basecalls = &calculateBasecalls($alignment);
		my $basecalllength=scalar @basecalls;
		
		
		# Save the raw base calls for future consensus calling
		# Write the base by base calls to a file....
		push(@{$data{uids}{$barcode}{raw}}, join("", @basecalls));
		

		# Store the individual basecalls in the readData hash
		for ( my $i = 0; $i < scalar(@basecalls); $i++ ) {
			my $basecall=$basecalls[$i];
			$basecall =~ tr/acgt/ACGT/;
			$data{basecalls}{$i}{raw}{$basecall}++;
			$data{basecalls}{$i}{raw}{depth}++;
		}
		

		
		
	}
	
	#print Dumper(%align_lengths);<STDIN>;
	
	return %data;
		
}



sub get_site_data_pe{
	
	### extracts read at a specific start point and stores data about these reads inot thhe has
	
	my($chrom,$start,$end,$sam,$AmpliconID,$offset)=@_;
	
	my %data;
	
	## GET ALL ALIGNMENTS THAT OVERLAP THIS that cover this positions
	#print STDERR "getting segments for $AmpliconID $chrom $start $end";<STDIN>;
	my $segment = $sam->segment($chrom,$start,$end);
	
	
	#my @alignments = $sam->features(-type=>'read_pair',-seq_id=>$chrom,-start=>$start,-end=>$end);
	
	#for my $alignment(@alignments){
		
	#	print Dumper($alignment);<STDIN>;
		
	#	my $read_name = $alignment->query->name();
	#	print "read_name=$read_name";<STDIN>;

		



	#}
	
	
	#my @alignments = $segment->features;
	## LIMIT TO ALIGNMENTS AT THIS EXACT START POSITION, and the paired end at the correct mate position	
    my @alignments = $segment->features(-filter => sub { my $alignment = shift;return ($alignment->start==$start) && ($alignment->mate_end==$end) });
	#my @alignments = $segment->features(-filter => sub { my $alignment = shift;return $alignment->start==$start });
	

	
	my $count_exact=scalar @alignments;
	my $count_offset=0;
	#if($offset){
	#	my @offset_alignments=$segment->features(-filter=>sub{ 
	#		my $alignment=shift;
	#		my $lower=$start-$offset;
	#		my $upper=$start+$offset;
	#		return $alignment if( ($alignment->start>=$lower) && ($alignment->start<=$upper) && ($alignment->start!=$start) );
	#	});
	#	$count_offset=scalar @offset_alignments;
	#	push(@alignments,@offset_alignments);
	#}

	### DEBUG, REMOVE for production
	#my $max=10000;
	#my $count=scalar @alignments;
	#$max=$count<$max ? $count : $max;
	#@alignments=@alignments[0...$max];	
	
	
	my $count_all=scalar @alignments;
	print STDERR "$count_all alignments for this site, $count_exact + $count_offset within $offset bases\n";

	### PROCESS EACH AlIGNMENT
	#my %align_lengths;
	for my $alignment(@alignments){
		
		
		#my $start=$alignment->start;
		#my $end=$alignment->end;
		#my $mstart=$alignment->mate_start;
		#my $mend=$alignment->mate_end;
		#print "$start $end $mstart $mend";<STDIN>;
		#next;
		
		#print "$inputSeqCount $chrom $chromStart";<STDIN>;
		$data{SitesSeqCount}++;
		

		
		my $barcode = '';
		my $bc_position = 0;
		my $read_name = $alignment->query->name();
		if ($read_name =~ /HaloplexHS/) { # We have a HaloplexHS read
			$read_name =~ m/HaloplexHS-([ACTG]{10})/;
			$barcode = $1;
		} else { # We have a SiMSenSeq read
			my $read_dna = $alignment->query->dna();
			($bc_position, $barcode) = &Debarcer::extractBarcodeQuick($read_dna);
		}
		next unless ( $bc_position == 0 );
		# Skip the barcode if it's in the mask hash loaded from the maskfile
		if ( exists $invalidBarcodes{$AmpliconID}->{$barcode} ) {
			# print "Skipping invalidBarcodes{$AmpliconID}->{$barcode}\n"; 
			next unless ( $opts{"justUIDdepths"} ); # Do not skip an invalid barcode if we're only generating the UID depth file.
		}
	
		# Count this instance of the barcode family
		$data{uids}{$barcode}{count}++;
	    
		### if not askig for basecalls, tehn no point going any futher
		next unless($opts{basecalls});
		
		my $ref_dna   = $alignment->dna; 
		my $query_dna = $alignment->query->dna; # query sequence bases
		#print "$query_dna\n";
		#print "$ref_dna\n";<STDIN>; 
		#$align_lengths{length($ref_dna)}++;
	
		# Determine the base call or insertion/deletion status wrt to the start of the alignment, i.e.
		# chromosomal position given in $chromStart.  Therefore $basecalls[0] is for $chromStart + 0
		my @basecalls = &calculateBasecalls($alignment);
		my $basecalllength=scalar @basecalls;
		
		
		# Save the raw base calls for future consensus calling
		# Write the base by base calls to a file....
		push(@{$data{uids}{$barcode}{raw}}, join("", @basecalls));
		

		# Store the individual basecalls in the readData hash
		for ( my $i = 0; $i < scalar(@basecalls); $i++ ) {
			my $basecall=$basecalls[$i];
			$basecall =~ tr/acgt/ACGT/;
			$data{basecalls}{$i}{raw}{$basecall}++;
			$data{basecalls}{$i}{raw}{depth}++;
		}
		

		
		
	}
	
	#print Dumper(%align_lengths);<STDIN>;
	
	return %data;
		
}


### this function will identify similar uids for a given amplicon,and either merge or mask
sub collapse_uids{
	
	my ($uids)=@_;
	return;
	
	my $levenshtien_binary = "$ENV{'BHOME'}/bin/levenshtien_stream";
	#print Dumper($uids);<STDIN>;
	
	my $uidcount=scalar keys %$uids;
	print "compare $uidcount uids\n";
	
	#(open my $TMP,">","uid.temp") || die "unable to open temp file";
	#print $TMP join("\n",sort keys %$uids);
	#close $TMP;
	
	#my @results=`cat uid.temp | $levenshtien_binary`;
	#my $count=scalar @results;
	#print "$count with editdist1";
	#<STDIN>;
	#exit;
	
	
	### capture in prefix and suffix hashes
	my %hash;
	for my $uid(sort {$$uids{$b}{count}<=>$$uids{$a}{count}} keys %$uids){
		
		
		
		my %testuids=map{$_=>1} grep{ ($$uids{$_}{count}/$$uids{$uid}{count}) < 0.2 }keys %$uids;
		my $testcount=scalar keys %testuids;
		print "$uid $$uids{$uid}{count} vs $testcount\n";
		
		
		
		my $uid_prefix=substr($uid, 0, 6);
		my $uid_suffix=substr($uid,length($uid)-6,6);
		push(@{$hash{prefix}{$uid_prefix}},$uid);
		push(@{$hash{suffix}{$uid_suffix}},$uid);
	}
	for my $xfix(qw/prefix suffix/){
		for my $subseq(keys %{$hash{$xfix}}){
			
			my @UIDs = @{$hash{$xfix}{$subseq}};
			my $nUIDs = scalar @UIDs;
			next if ( $nUIDs == 1 );
			
			print "assessng $nUIDs with $xfix $subseq\n";
			
			my $UIDstring = join("\n", @UIDs);
			my @closeUIDs = `echo "$UIDstring" | $levenshtien_binary`;
			next unless ( scalar(@closeUIDs) );
			chomp @closeUIDs;
			#print "$subseq\n";
			#print Dumper(@UIDs);<STDIN>;
			#print Dumper(@closeUIDs);<STDIN>;
		
			for (@closeUIDs)  {
				next unless (/^Edit/);
				my ( $descriptor, $uid1, $uid2 ) = split("\t", $_);
				my $d1 = $$uids{$uid1}{count};
				my $d2 = $$uids{$uid2}{count};
				my ($uid,$identity,$uid_depth,$identity_depth,$ratio)=$d1>$d2 ? ($uid2,$uid1,$d2,$d2/$d1) : ($uid1,$uid2,$d1,$d1/$d2);
				if ( ($ratio < 0.2) && ($uid_depth > 5) ) { # This is the masking heuristic
					$$uids{$uid}{collapse}{mask}=1;
					$$uids{$uid}{collapse}{identity}=$identity;
					$$uids{$uid}{collapse}{ratio}=$ratio;
					$$uids{$uid}{collapse}{identity_depth}=$identity_depth;
				}
			}
		}
	}
}







## function to make a call table at specific depth thresholds
sub get_position_calls{
	
	
	my ($id,$data,$thresholds)=@_;
	### thresholds is a comma separated string listing thresholds
	my @thresholds=split /,/,$thresholds;
	my $threshold=shift @thresholds;
	
	print STDERR "assessing positions calls as threshold $threshold\n";
	
	### $data is a reference to the site data
	my @SNVtypes=qw/A C G T D I N/; 
	print STDERR "calculating position table for $id\n";
		
	my %table;
	

	foreach my $position ( sort { $a <=> $b } keys %{$$data{basecalls}}) {
		my ($probableBase, $n, $rawDataString) = ('', 0, '');
		
		print "position=$position\n";
		#print Dumper($$data{basecalls}{$position});<STDIN>;
		
		### for this position, assess the raw data
		
		
		### I think this can be omitted, it is not saved at all.  ie, probable base information.  and it is recalculated in the main function
		### check on this to ensure you've not missed any logic
		my %rawDataHash = ();
		my $depth = 0;  #Depth
		for my $snv ( @SNVtypes ) {
			# Save the base identity if the count is higher than any previous observed count
			my $snvdepth=$$data{basecounts}{$position}{$snv} || 0;
			### probable base is the snv with the hightest depth
			if ( $snvdepth > $n ) {
				$probableBase = $snv;
				$n = $snvdepth;
			}
				
			$rawDataString .= "\t" . $snvdepth;
			$rawDataHash{$snv} = $snvdepth;
			$depth += $snvdepth;
		}

		# Do not report this site if it's a low coverage position
		next if ( $depth < $threshold );
		
		# Calculate the real genome position
		my ($chrom, $chromStart) = split(/:/, $id);
		my $genomePosition = $chromStart + $position;

			# This section will restrict reporting of amplicon positions to those within
			# the target window specified in the *amplicons.txt file
#			if ( $opts{"justTargets"} ) {
#				# If a target window exists, skip $genomePosition unless it falls within start and end
#				if ( exists $ampliconInfo{$ampliconName}{"TargetWindow"} ) {
#					unless ( $genomePosition >= $ampliconInfo{$ampliconName}{"TargetWindow"}{"start"} & $genomePosition <= 	$ampliconInfo{$ampliconName}{"TargetWindow"}{"end"} ) {
#						# print STDERR "Skipping $ampliconName $amp Position $position $genomePosition\n";
#						next;
#					}
#				}
#				# To gen here, no target window exists, so report everything
#			}
		
			# Save accumulated data in the output table
	
		foreach my $base ( @SNVtypes ) {
			$table{$genomePosition}{"raw"}{$base} += $rawDataHash{$base};
			$table{$genomePosition}{"consensus"}{$base} += $$data{conspos}{$position}{$base};
		}
	}
	return %table;
	#print "ready to dump table";<STDIN>;
	#print Dumper(%table);<STDIN>;
	
	
}






## collapse all the reads for a given uuid to a consensus sequence
sub generate_consensus_data{
	my($id,$FH,$hashref,$familysizes)=@_;
	
	my @levels=split /,/,$familysizes;
	my %uids=%{$$hashref{uids}};
	my $uidcount=scalar keys %uids;
	#print STDERR "generating consensus data for $id : $uidcount uids\n";
	

	
	my ($AmpliconCount, $AmpliconCoverage) = (0, 0);
	### get amplicon coverage at this depth
	for my $uid ( sort {$uids{$b}{count}<=>$uids{$a}{count}} keys %uids) {
		my $uidcount=$uids{$uid}{count};
		
		my @uidlevels=grep{$uidcount>=$_} @levels;
		#print "$uid $uidcount " . join(",",@uidlevels);<STDIN>;
		
		$AmpliconCount++;
		$AmpliconCoverage += $uidcount;
		

		### store the consensus back into the data reference
		my @raw_reads=@{$uids{$uid}{raw}};
		my $raw_read_count=scalar @raw_reads;
		#print "$id $uid $raw_read_count";<STDIN>;
		
		## limit consensus to 100 raw_reads?  why use all of them, sometimes > 1000.  ?????
		#@raw_reads=@raw_reads[0..99] if($raw_read_count>100);
		
		### get the consenus sequecne, with a minium depth. if not the miniumum depth then no consensus is regturned
		my $consensus=&generateConsensus(@raw_reads,1) || "N/A";   ### second argument used to be miniumum depth, now set to 1, do all depths
		
		
		#print Dumper(@raw_reads);<STDIN>;
		#print "consensus $consensus";<STDIN>;
		
		
		### store the consensus sequence back into the uids sturcture reference here
		$$hashref{uids}{$uid}{consensus}=$consensus;   
		
		# $familyData{$amp}->{$barcode}{"consensus"} = &generatePhyloConsensus(@{$familyData{$amp}->{$barcode}{"raw"}}, $consensusDepth);  # Still in progress.
		# Write the UID depth information to a file
		printf $FH ("%s\t%s\t%s\t%s\n", $id, $uid, $uidcount, $consensus );
		
		
		
		
		#### now organize the bases at this position into a hash structure
		my @bases = split(//, $consensus);
		for (my $i = 0; $i < scalar(@bases); $i++) {
			my $base=$bases[$i];
			$base =~ tr/acgt/ACGT/;
			
			for my $uidlevel(@uidlevels){
				$$hashref{basecalls}{$i}{consensus}{$uidlevel}{$base}++;
				$$hashref{basecalls}{$i}{consensus}{$uidlevel}{depth}++;
			}	
		}
	}
	
	
	
	print STDERR "$id\tcount|coverage\t$AmpliconCount\t$AmpliconCoverage\n";
}

sub print_uid_depths{
	my ($FH,$id,$hashref)=@_;
	for my $uid(sort { $$hashref{$b}{count}<=>$$hashref{$a}{count} } keys %$hashref){
		printf $FH ("%s\t%s\t%s\n", $id, $uid, $$hashref{$uid}{"count"} );
	}
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
	#print "Ref:  $ref\n      $matches\nRead: $query\n\n";
	#<STDIN>;
	# Since query is longer than ref due to adapters, etc.,
	# trim the sequences
	$ref =~ /^(-+).+?(-+)$/;
	
	($ref, $query) = ( substr($ref, length($1), length($ref) - (length($1 . $2)) ), substr($query, length($1), length($query) - (length($1 . $2)) ) );
	($matches) = ( substr($matches, length($1), length($matches) - (length($1 . $2)) ) );
	#print "Mismatch:\n" if ( $debug & ($ref ne $query) );
	#print "Ref:  $ref\n      $matches\nRead: $query\n"; # if ( $debug );
	#<STDIN>;
	
	# <STDIN>;  # Pause for keypress
	
	# Index the sequence by genomic position
	my $genomicOffset = 0;
	my $inInsertion = 0;
	
	my @Ar = split(//, $ref);    #print Dumper(@Ar);<STDIN>;
	my @Aq = split(//, $query);  #print Dumper(@Aq);<STDIN>;
		
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
	
	#print join("", "CIGAR:", @basecalls, "*\n\n");  #if ( $debug );
	#<STDIN>;
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

sub calculateAmpliconReportingThreshold_rev {

=pod

=head2 calculateAmpliconReportingThreshold

This function identifies a minimal threshold for each amplicon,
below which positions aren't reported in the cons<depth>.txt files

=cut


	### appears to be finding the maximum depth across ALL positions, and setting depthcut to 10% of that value
	my %href = @_;  # This is %readData
	# readData format is
	# $readData{$position}{$basecall};
	my $depthCut;
	
	# find the maximum depth for each amplicon ( usually the first site )
	foreach my $position ( sort {$a<=>$b} keys %href ) {
		#print Dumper($href{$position});<STDIN>;
		my $depthHere = 0;
		$depthHere += $href{$position}{$_} foreach ( keys %{$href{$position}} );
		$depthCut = $depthHere if ( $depthHere > $depthCut );
		#print STDERR "$position $depthHere $depthCut";<STDIN>;

	}
	
	# Adjust depth cuts downward to a percentage of max
	$depthCut = int( $depthCut * 0.1 ); 
	
	return $depthCut;
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

sub validate_options{
	my ($opts)=@_;
	usage("Help requested.") if($opts{help});
	
	if(! $$opts{configfile} || ! -e $$opts{configfile}){
		usage("Configuration file not provided or not found.");
	}
	if(! $$opts{bam} || ! -e $$opts{bam}){
		usage("Bam file not provided or not found.");
	}
	
	if(! $$opts{sampleID}){
		usage("sample ID not provided.");
	}

	if(! $$opts{output_folder} || ! -d $$opts{output_folder}){
		usage("Output folder not provided or not found.");
	}
	
	if(! $opts{refgenome} || -d $opts{refgenome}){
		usage("Reference Genome not provided in the configuration file, or file not found");
	}
	$$opts{consDepth}=3 unless($$opts{consDepth});
	$$opts{plexity}=1 unless($$opts{plexity});
	$$opts{offset}=0 unless($$opts{offset});
}

sub usage{
	print "\ndebarc.pl [options]\n";
	print "Options are as follows:\n";
	print "\t--config String/filename. \n";
	print "\t--bam String/filename.  Aligned sequence data with uid information embedded in the header lines.\n";
	print "\t--sampleID String. Name used to prefix output files.\n";
	print "\t--consDepth Integer. The mininum family size required to calculates consensus information. Default depth is 3\n";
	print "\t--plexity Integer. The number of sites to assess.  Will emit to sites file if it does not exist, or process this many sites. Default is 1.\n";
	print "\t--sites String/File.  A file to write sites (if file doesn't exist), or read sites that will be assessed.\n";
	print "\t--site String.  A single site to assess. Cannot be used in combination with -sites\n";
	print "\t--UIDdepths.  Emit the UIDdepths to a file.\n";
	print "\t--basecalls.  Calculate consensus and base counts at each position.\n";
	print "\t--out.  String/directory. Where to save analysis. Subfolders : tables, : will be created in this folder\n";
	print "\t--help displays this usage message.\n";

	die "\n@_\n\n";
}



