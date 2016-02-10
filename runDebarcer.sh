

# Generate the files used by R analysis script to summarize in a table.
#

if [[ $# -ne 2 ]]; then
	echo "
Need to specify a filename and samplename as arguments.
	
	Usage: runDebarcer.sh <infile.fastq.gz> <SampleName>
	";
	exit 1;
fi

# . /etc/profile
# module load samtools/0.1.19;

# Uncomment to just regenerate the graphics
# . generateGraphicalReports.sh $BHOME $2; exit;  # $2 is samplename

SVNID="$Id: generatePositionalCompositionFiles.sh 387 2016-02-09 22:08:35Z pkrzyzanowski $"
STARTTIME=$(date +"%Y%m%d-%H%M%S")
MAINLOG="DeBarcEr."$STARTTIME".log"
echo "[Debarcer `date`] : Logfile timestamp: $STARTTIME" > $MAINLOG
echo "Running Debarcer version $SVNID" >> $MAINLOG
echo "Running in: `pwd`" >> $MAINLOG

mkdir -p sge # For child process log files

# FASTQGZ="/u/pkrzyzanowski/projects/EAC/molecular_barcoding/SaferSeqTests/14July02/fastqs/Sample_01.R1.fastq.gz"
# FASTQGZ="/u/pkrzyzanowski/projects/EAC/molecular_barcoding/data/genomics.med.tufts.edu/140721-070_Jennifer_Jackson/sequence_data_illumina/Unaligned/Sample_6.R1.fastq.gz"
# SAMPLENAME="Sample01"
FASTQGZ=$1
SAMPLENAME=$2

# Defaults, over-rideable by a local debarcer.config file.
typeset -A config
config=(
	[plexity]=5
	)

if [ -e debarcer.conf ]; then
	while read line
	do
		if echo $line | grep -F = &>/dev/null
		then
			varname=$(echo "$line" | cut -d '=' -f 1)
			config[$varname]=$(echo "$line" | cut -d '=' -f 2-)
			echo "Resetting $varname to "${config[$varname]} >> $MAINLOG;
		fi
	done < debarcer.conf
fi

PLEXITY=${config[plexity]}
echo "PLEXITY=$PLEXITY" >> $MAINLOG

# Test added so that $BHOME can be set as an environment variable to allow for use of exported copies
# in analysis directories
# This is set by the module file until ./utils when the module is loaded
if [ ! "$BHOME" ]; then
	BHOME="/u/pkrzyzanowski/projects/EAC/molecular_barcoding/svn/trunk/debarcer"
fi
echo "[Debarcer `date`] Running workflow from $BHOME" >> $MAINLOG

# Optional setting
AMPLICON_TABLE="$BHOME/amplicon_tables/all_amplicons.txt";  # This is default

# Some bwa mem files for future development
#
if [ ! -e $SAMPLENAME.$FASTQGZ.sorted.bam.touch ]; then
	$BHOME/tools/runBWA.sh $FASTQGZ $SAMPLENAME
	touch $SAMPLENAME.$FASTQGZ.sorted.bam.touch
fi
echo "[Debarcer `date`] Raw reads mapped by bwa: `samtools view $SAMPLENAME.$FASTQGZ.sorted.bam | wc -l`" >> $MAINLOG
# samtools view $SAMPLENAME.$FASTQGZ.sorted.bam | cut -f 3,4 | perl $BHOME/tools/uniqCount.pl | tail -n 20 >> $MAINLOG # List top 20 amplicons, for testing

# EXPERIMENTAL SECTION, STILL IN DEVELOPMENT
# time $BHOME/generateDownsamplingEstimates.sh $SAMPLENAME.$FASTQGZ.sorted.bam $SAMPLENAME
#
# End test section


# This section is not needed anymore since the raw counts are 
# always output in the generateConsensusFromBAM.pl tables
#
# echo "[Debarcer `date`] Running $SAMPLENAME Raw" >> $MAINLOG
# if [ ! -e $SAMPLENAME.rawPositionalComposition.txt.touch ]; then
# time perl $BHOME/classifyReadsToAmplicons.pl --bamfile=$SAMPLENAME.$FASTQGZ.sorted.bam --dump --ampTable=$AMPLICON_TABLE 2> >(tee -a $MAINLOG >&2) | 
    # perl $BHOME/reportCommonPositionBases.pl --raw > $SAMPLENAME.rawPositionalComposition.txt 2> >(tee -a $MAINLOG >&2)
	# touch $SAMPLENAME.rawPositionalComposition.txt.touch
# fi

# This was supposed to generate a consensus sequences file
#    
# echo "[Debarcer `date`] Running $SAMPLENAME Consensus Sequences" >> $MAINLOG
# if [ ! -e $SAMPLENAME.consensusSequences.txt.touch ]; then
# time perl $BHOME/classifyReadsToAmplicons.pl --bamfile=$SAMPLENAME.$FASTQGZ.sorted.bam --consensus --strictCons --ampTable=$AMPLICON_TABLE |  # Set strict setting always
	# cat | gzip > $SAMPLENAME.consensusSequences.txt.gz
	# touch $SAMPLENAME.consensusSequences.txt.touch
# fi


echo "[Debarcer `date`] Generating UID depth file for $SAMPLENAME" >> $MAINLOG
rm -f $SAMPLENAME.barcode_mask # Remove the mask file prior to identifying masked barcodes
time perl $BHOME/generateConsensusFromBAM.pl --bam=$SAMPLENAME.$FASTQGZ.sorted.bam --sampleID=$SAMPLENAME --plexity=$PLEXITY --justUIDdepth 2> >(tee -a $MAINLOG >&2)

# Creation of a barcode masking script goes here.
echo "[Debarcer `date`] Creating barcode mask file" >> $MAINLOG
gunzip -c $SAMPLENAME.UIDdepths.txt.gz | perl $BHOME/tools/identifyMaskableBarcodes.pl > $SAMPLENAME.barcode_mask
#
# Comment: subsequent runs of generateConsensusFromBAM.pl will use the mask file
# to regenerate the UID.depths file, without the masked barcodes.


echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=1" >> $MAINLOG
# Arguments: --sampleID; --consDepth; --plexity ... others.
qsub -N DbC1$SAMPLENAME -l h_vmem=16G -e DbC1$SAMPLENAME.log -o DbC1$SAMPLENAME.log -cwd -b y "module load debarcer/trunk; time perl $BHOME/generateConsensusFromBAM.pl --bam=$SAMPLENAME.$FASTQGZ.sorted.bam --sampleID=$SAMPLENAME --consDepth=1 --plexity=$PLEXITY > $SAMPLENAME.bamPositionComposition.cons1.txt"
	
echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=3" >> $MAINLOG
# Arguments: --sampleID; --consDepth; --plexity ... others.
qsub -N DbC3$SAMPLENAME -l h_vmem=16G -e DbC3$SAMPLENAME.log -o DbC3$SAMPLENAME.log -cwd -b y "module load debarcer/trunk; time perl $BHOME/generateConsensusFromBAM.pl --bam=$SAMPLENAME.$FASTQGZ.sorted.bam --sampleID=$SAMPLENAME --consDepth=3 --plexity=$PLEXITY > $SAMPLENAME.bamPositionComposition.cons3.txt"
	
echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=10" >> $MAINLOG
# Arguments: --sampleID; --consDepth; --plexity ... others.
qsub -N DbC10$SAMPLENAME -l h_vmem=16G -e DbC10$SAMPLENAME.log -o DbC10$SAMPLENAME.log -cwd -b y "module load debarcer/trunk; time perl $BHOME/generateConsensusFromBAM.pl --bam=$SAMPLENAME.$FASTQGZ.sorted.bam --sampleID=$SAMPLENAME --consDepth=10  --plexity=$PLEXITY > $SAMPLENAME.bamPositionComposition.cons10.txt"
	
echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=20" >> $MAINLOG
# Arguments: --sampleID; --consDepth; --plexity ... others.
qsub -N DbC20$SAMPLENAME -l h_vmem=16G -e DbC20$SAMPLENAME.log -o DbC20$SAMPLENAME.log -cwd -b y "module load debarcer/trunk; time perl $BHOME/generateConsensusFromBAM.pl --bam=$SAMPLENAME.$FASTQGZ.sorted.bam --sampleID=$SAMPLENAME --consDepth=20  --plexity=$PLEXITY > $SAMPLENAME.bamPositionComposition.cons20.txt"
	
echo "[Debarcer `date`] BAM Consensus for $SAMPLENAME depth	=30" >> $MAINLOG
# Arguments: --sampleID; --consDepth; --plexity ... others.
qsub -N DbC30$SAMPLENAME -l h_vmem=16G -e DbC30$SAMPLENAME.log -o DbC30$SAMPLENAME.log -cwd -b y "module load debarcer/trunk; time perl $BHOME/generateConsensusFromBAM.pl --bam=$SAMPLENAME.$FASTQGZ.sorted.bam --sampleID=$SAMPLENAME --consDepth=30  --plexity=$PLEXITY > $SAMPLENAME.bamPositionComposition.cons30.txt"

qsub -N "AggregateDebarcers" -hold_jid DbC1$SAMPLENAME,DbC3$SAMPLENAME,DbC10$SAMPLENAME,DbC20$SAMPLENAME,DbC30$SAMPLENAME -sync y -cwd -b y -e ./sge -o ./sge "sleep 1"

cat DbC*.log >> $MAINLOG  # Aggregate the log files generated above.
rm DbC*.log


# generate nucleotide context error tables
# Still in development
#
# perl $BHOME/tools/trinucleotideErrors.pl --infile=$SAMPLENAME.rawPositionalComposition.txt > $SAMPLENAME.trinucErrors.raw.txt
# perl $BHOME/tools/trinucleotideErrors.pl --infile=$SAMPLENAME.consensusPositionalComposition.bc10.txt > $SAMPLENAME.trinucErrors.cons10.txt
# perl $BHOME/tools/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.rawPositionalComposition.txt > $SAMPLENAME.trinucErrors.raw.perSite.txt
# perl $BHOME/tools/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.consensusPositionalComposition.bc10.txt > $SAMPLENAME.trinucErrors.cons10.perSite.txt
# perl $BHOME/tools/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.rawPositionalComposition.txt --pentanucleotides > $SAMPLENAME.pentanucErrors.raw.perSite.txt
# perl $BHOME/tools/trinucleotideErrors-perPosition.pl --infile=$SAMPLENAME.consensusPositionalComposition.bc10.txt --pentanucleotides > $SAMPLENAME.pentanucErrors.cons10.perSite.txt
 

echo "[Debarcer `date`] Running $SAMPLENAME Barcode Distribution Report" >> $MAINLOG # Approx 1 minute runtime
if [ ! -e $SAMPLENAME.barcodeComposition.txt.touch ]; then
time gunzip -c $SAMPLENAME.UIDdepths.txt.gz |
	perl $BHOME/reportBarcodeComposition.pl > $SAMPLENAME.barcodeComposition.txt
	touch $SAMPLENAME.barcodeComposition.txt.touch
fi

# Generate the graphics.
. generateGraphicalReports.sh $BHOME $2;


# Generate summary statistics file
cat $MAINLOG | perl $BHOME/tools/summarizeAmpliconYields.pl --sampleID=$SAMPLENAME > $SAMPLENAME.SummaryStatistics.txt

gunzip -c $SAMPLENAME.UIDdepths.txt.gz | perl $BHOME/tools/summarizeAmpliconConsensusDepths.pl --sampleID=$SAMPLENAME --depths=1,3,10,20,30,100 > $SAMPLENAME.consensusStatistics.txt

# Package up results files
tar cvfz $SAMPLENAME.DebarcerResultsPackage.$STARTTIME.tar.gz `find . -name \"*pdf\" -o -name \"*.log\" -o -name \"*Statistics*\" -o -name \"*PositionComposition*\"`


# SVNID: $Id: generatePositionalCompositionFiles.sh 387 2016-02-09 22:08:35Z pkrzyzanowski $
