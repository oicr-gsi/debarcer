

# Generate the downsampling tables using the bam file with aligned amplicons.
#

if [[ $# -ne 2 ]]; then
	echo "
Need to specify a filename and samplename as arguments.
	
	Usage: generateDownsamplingEstimates.sh <infile.bam> <SampleName> <DownsampleLevel>
	";
	exit 1;
fi

mkdir -p downsampling/sgelog;

BAMFILE=$1
SAMPLENAME=$2
DS=1  # Downsampling decile for samtools.  i.e. 1 = 10%

for DS in `seq 9`; do
	for rep in `seq 10`; do
	qsub -N "bcDownsample" -l h_vmem=16G -e ./downsampling/sgelog -o ./downsampling/sgelog -cwd -b y "module load samtools/0.1.19; export BHOME=$BHOME; time perl $BHOME/generateConsensusFromBAM.pl --bam=$BAMFILE --downsample=$rep.$DS | perl ../BOARAT/tools/summarizeConsensusCounts.pl --downsample=0.$DS > ./downsampling/$SAMPLENAME.downsampling.0.$DS.replicate$rep.txt";
	done;	
done;	


# To plot, aggregate all results and in R:
#
# ds <- read.delim(file="downsampling.txt", header=F)
# ds <- aggregate(V3 ~ V1 + V2 + V4, data = ds, mean)
# qplot(V4, V3, data=ds, colour=factor(V2)) + geom_jitter()
# or
# qplot(V4, V3, data=ds[ds$V2 == 10,], colour=factor(V1), geom="line")
# qplot(V4, V3, data=ds[ds$V2 == 3,], colour=factor(V1), geom="smooth") # without aggregate
# qplot(V4, V3, data=ds[ds$V1 == "chr19:11132395",], colour=factor(V2), geom="smooth") # without aggregate








