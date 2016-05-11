

# Generate the downsampling tables using the *UIDdepths.txt.gz file.

SAMPLENAME=$1
Rscript $BHOME/src/quickDownsamplingEstimates.R ./tables/$SAMPLENAME.UIDdepths.txt.gz