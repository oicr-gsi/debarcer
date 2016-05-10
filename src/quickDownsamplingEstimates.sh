

# Generate the downsampling tables using the *UIDdepths.txt.gz file.

mkdir -p downsampling;

Rscript $BHOME/tools/quickDownsamplingEstimates.R


# To plot, aggregate all results and in R:
#
# ds <- read.delim(file="downsampling.txt", header=F)
# ds <- aggregate(V3 ~ V1 + V2 + V4, data = ds, mean)
# qplot(V4, V3, data=ds, colour=factor(V2)) + geom_jitter()
# or
# qplot(V4, V3, data=ds[ds$V2 == 10,], colour=factor(V1), geom="line")
# qplot(V4, V3, data=ds[ds$V2 == 3,], colour=factor(V1), geom="smooth") # without aggregate
# qplot(V4, V3, data=ds[ds$V1 == "chr19:11132395",], colour=factor(V2), geom="smooth") # without aggregate








