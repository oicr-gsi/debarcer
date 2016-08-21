#!/bin/bash

# Debarcer
#
# This is the main run script for the Debarcer pipeline
# 
# Debarcer aligns, demultiplexes and analyzes barcoded
# NGS sequencing data in the form of fastq files, and
# will generate reports in the directory from which
# it is run.
#
# For more information, see the README files in the
# docs directory.
#
# Author: Paul Krzyzanowski
# Email:  paul.krzyzanowski@oicr.on.ca
# (c) 2014-2016

FASTQGZ=''
SAMPLENAME=''
LAUNCHDIR=$PWD
OUTPUTDIR=$LAUNCHDIR
DBROOT=$BHOME
DBSRC=$DBROOT"/src/"
USE_SGE=false

usage()
{
	echo "
Need to specify a run mode, filename, samplename, and output directory as arguments.

	Usage: runDebarcer.sh [-u|-g|-r] [Options] -f <infile.fastq.gz> -n <SampleName> -o <output dir>
";
}

# Check for any options being passed
if [[ ! $1 ]]; then
	usage;
	echo "Error! No options were passed!";
	exit;
fi

while getopts ":gruf:n:o:" opt; do
	case $opt in
		u)
			usage;
			exit ;1;
			;;
		r)
			echo "Running debarcer..."; >&2
			;;
		g)
			ONLYGRAPHICS=1
			;;
		f)
			FASTQGZ=$OPTARG;
			FASTQBASENAME=${FASTQGZ##*/}
			;;
		n)
			SAMPLENAME=$OPTARG;
			;;
		o)
			OUTPUTDIR=$OPTARG;
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1;
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			echo "Try runDebarcer.sh -u for usage";
			exit;
			;;
	esac
done

mkdir -p $OUTPUTDIR
cd $OUTPUTDIR

SAMPLEPREFIX=$SAMPLENAME"."$FASTQBASENAME

# Generate graphics only
if [[ $ONLYGRAPHICS ]]; then
	echo "Debarcer: Generating graphics only." >&2
	. generateGraphicalReports.sh $DBROOT $SAMPLENAME; 
	exit;
fi


VERSIONID="0.4.0"
STARTTIME=$(date +"%Y%m%d-%H%M%S")
MAINLOG="DeBarcEr."$STARTTIME".log"
PERLVERSION=`perl -e 'print $^V';`
echo "[Debarcer `date`] : Logfile timestamp: $STARTTIME" > $MAINLOG
echo "Running Debarcer version $VERSIONID" >> $MAINLOG
echo "Running perl version $PERLVERSION" >> $MAINLOG
echo "---" >> $MAINLOG
echo "git info:" >> $MAINLOG
cat $DBROOT/.git/HEAD >> $MAINLOG
cat $DBROOT/.git/`cut -d" " -f 2 $DBROOT/.git/HEAD` >> $MAINLOG
echo "---" >> $MAINLOG

echo "Running in: `pwd`" >> $MAINLOG

echo 'Sourcing ~/.debarcer'
source ~/.debarcer

mkdir -p sge # For child process log files
mkdir -p tables
mkdir -p figures


# Test added so that $BHOME can be set as an environment variable to allow for use of exported copies
# in analysis directories
# This is set by the module file from ./utils/modulefiles/debarcer/ when the module is loaded
if [ ! "$BHOME" ]; then
	echo 'You need to set the $BHOME environment variable to the debarcer location.';
	echo 'If modules are being used, do "module load debarcer"';
	exit 1;	
fi
echo "[Debarcer `date`] Running workflow from $BHOME" >> $MAINLOG

CONFIG_FILE="$DBROOT/config/debarcer.conf"
echo "Using debarcer config file:" >> $MAINLOG
cat $CONFIG_FILE >> $MAINLOG;
echo "End of debarcer config file." >> $MAINLOG

# Some bwa mem files for future development
#
if [ ! -e $SAMPLEPREFIX.sorted.bam.touch ]; then
	echo "[Debarcer `date`] Running runBWA.sh" >> $MAINLOG;
	$DBSRC/runBWA.sh $FASTQGZ $SAMPLENAME;
	touch $SAMPLEPREFIX.sorted.bam.touch;
fi
echo "[Debarcer `date`] Raw reads mapped by bwa: `$SAMTOOLSROOT/bin/samtools view $SAMPLEPREFIX.sorted.bam | wc -l`" >> $MAINLOG


echo "[Debarcer `date`] Generating UID depth file for $SAMPLENAME" >> $MAINLOG
rm -f ./tables/$SAMPLENAME.barcode_mask # Remove the mask file prior to identifying masked barcodes
time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --config=$CONFIG_FILE --justUIDdepth 2> >(tee -a $MAINLOG >&2)

# Creation of a barcode masking script goes here.
echo "[Debarcer `date`] Creating barcode mask file" >> $MAINLOG
gunzip -c ./tables/$SAMPLENAME.UIDdepths.txt.gz | perl $DBSRC/identifyMaskableBarcodes.pl > ./tables/$SAMPLENAME.barcode_mask
#
# Comment: subsequent runs of generateConsensusFromBAM.pl will use the mask file
# to regenerate the UID.depths file, without the masked barcodes.


# Generate downsampling estimator plots
time $BHOME/src/quickDownsamplingEstimates.sh $SAMPLENAME
rm Rplots.pdf # Remove empty plot just made by ggplot


## FIXME Since the UID depths file is available the composition can be analyzed here using src/_calculateBarcodeBiases.pl


if [ $USE_SGE = true ]; then

	echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=1" >> $MAINLOG
	# Arguments: --sampleID; --consDepth; --plexity ... others.
	qsub -N DbC1$SAMPLENAME -l h_vmem=16G -e DbC1$SAMPLENAME.log -o DbC1$SAMPLENAME.log -cwd -b y "module load debarcer; time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=1 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons1.txt"
		
	echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=3" >> $MAINLOG
	# Arguments: --sampleID; --consDepth; --plexity ... others.
	qsub -N DbC3$SAMPLENAME -l h_vmem=16G -e DbC3$SAMPLENAME.log -o DbC3$SAMPLENAME.log -cwd -b y "module load debarcer; time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=3 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons3.txt"
		
	echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=10" >> $MAINLOG
	# Arguments: --sampleID; --consDepth; --plexity ... others.
	qsub -N DbC10$SAMPLENAME -l h_vmem=16G -e DbC10$SAMPLENAME.log -o DbC10$SAMPLENAME.log -cwd -b y "module load debarcer; time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=10 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons10.txt"
		
	echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=20" >> $MAINLOG
	# Arguments: --sampleID; --consDepth; --plexity ... others.
	qsub -N DbC20$SAMPLENAME -l h_vmem=16G -e DbC20$SAMPLENAME.log -o DbC20$SAMPLENAME.log -cwd -b y "module load debarcer; time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=20 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons20.txt"
		
	echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=30" >> $MAINLOG
	# Arguments: --sampleID; --consDepth; --plexity ... others.
	qsub -N DbC30$SAMPLENAME -l h_vmem=16G -e DbC30$SAMPLENAME.log -o DbC30$SAMPLENAME.log -cwd -b y "module load debarcer; time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=30 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons30.txt"

	qsub -N "AggregateDebarcers" -hold_jid DbC1$SAMPLENAME,DbC3$SAMPLENAME,DbC10$SAMPLENAME,DbC20$SAMPLENAME,DbC30$SAMPLENAME -sync y -cwd -b y -e ./sge -o ./sge "sleep 1"

	cat DbC*.log >> $MAINLOG  # Aggregate the log files generated above.
	rm DbC*.log

else

	# Do the commands serially on the local node
	time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=1 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons1.txt
	time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=3 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons3.txt
	time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=10 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons10.txt
	time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=20 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons20.txt
	time perl $DBSRC/generateConsensusFromBAM.pl --bam=$SAMPLEPREFIX.sorted.bam --sampleID=$SAMPLENAME --consDepth=30 --config=$CONFIG_FILE > ./tables/$SAMPLENAME.bamPositionComposition.cons30.txt

fi

# FIXME generate nucleotide context error tables
# Still in development
#
# perl $DBSRC/trinucleotideErrors.pl --infile=$SAMPLENAME.rawPositionalComposition.txt > $SAMPLENAME.trinucErrors.raw.txt
# perl $DBSRC/trinucleotideErrors.pl --infile=$SAMPLENAME.consensusPositionalComposition.bc10.txt > $SAMPLENAME.trinucErrors.cons10.txt
# perl $DBSRC/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.rawPositionalComposition.txt > $SAMPLENAME.trinucErrors.raw.perSite.txt
# perl $DBSRC/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.consensusPositionalComposition.bc10.txt > $SAMPLENAME.trinucErrors.cons10.perSite.txt
# perl $DBSRC/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.rawPositionalComposition.txt --pentanucleotides > $SAMPLENAME.pentanucErrors.raw.perSite.txt
# perl $DBSRC/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.consensusPositionalComposition.bc10.txt --pentanucleotides > $SAMPLENAME.pentanucErrors.cons10.perSite.txt
 

echo "[Debarcer `date`] Running $SAMPLENAME Barcode Distribution Report" >> $MAINLOG # Approx 1 minute runtime
if [ ! -e ./tables/$SAMPLENAME.barcodeComposition.txt.touch ]; then
time gunzip -c ./tables/$SAMPLENAME.UIDdepths.txt.gz |
	perl $DBSRC/reportBarcodeComposition.pl > ./tables/$SAMPLENAME.barcodeComposition.txt
	touch ./tables/SAMPLENAME.barcodeComposition.txt.touch
fi

# Generate the graphics.
. generateGraphicalReports.sh $DBROOT $SAMPLENAME;


# Generate summary statistics files
# These should stay in the root results directory
cat $MAINLOG | perl $DBSRC/summarizeAmpliconYields.pl \
	--config=$CONFIG_FILE  \
	--sampleID=$SAMPLENAME > $SAMPLENAME.SummaryStatistics.txt
gunzip -c ./tables/$SAMPLENAME.UIDdepths.txt.gz | 
	perl $DBSRC/summarizeAmpliconConsensusDepths.pl --sampleID=$SAMPLENAME --depths=1,3,10,20,30,100 > $SAMPLENAME.consensusStatistics.txt

cd $LAUNCHDIR
exit
