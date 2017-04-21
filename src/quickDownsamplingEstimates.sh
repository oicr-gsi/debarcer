

# Generate the downsampling tables using the *UIDdepths.txt.gz file.
module load R/3.3.0

SAMPLENAME=$1
Rscript $BHOME/src/quickDownsamplingEstimates.R ./tables/$SAMPLENAME.UIDdepths.txt.gz