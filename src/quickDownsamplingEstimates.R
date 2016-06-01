


UIDdepthFile <- as.character(commandArgs(trailingOnly = TRUE)[1])
outputdir <- "./figures/"

paste("Running bam_analysis.ErrorBarplots.txt")

library(ggplot2)

# Load data and check contents
UIDs <- read.delim(file = UIDdepthFile, header = F)
colnames(UIDs) <- c("Amplicon", "Barcode", "Count")
# str(UIDs)

# Determine how many amplicons of each coordinate are present
allAmplicons <- as.character(unique(UIDs$Amplicon))
ampliconCounts <- sort(table(UIDs$Amplicon), decreasing = TRUE)

a <- names(ampliconCounts)[1]  # Select most common amplicon

UIDs.a <- UIDs[UIDs$Amplicon == a,]

results <- data.frame(depth = 0, barcodes = 0)

# Define downsampling cutoffs
actualDepth <- sum(UIDs$Count)
cutoffs <- round(seq(0, 1, 0.1) * actualDepth)
print("Using depth cutoffs:")
cutoffs

for ( d in cutoffs ) {
	sample_depth <- d

	estimateBarcodes <- function() {
		UIDs.a$prob <- UIDs.a$Count / sum(UIDs.a$Count)
		UIDs.a$BarcodesAsChar <- as.character(UIDs.a$Barcode)

		sample_barcodes <- sample(UIDs.a$BarcodesAsChar, size = sample_depth, replace = T, prob = UIDs.a$prob )
		sample_barcodes.length <- length(table(sample_barcodes))
		return(sample_barcodes.length)
	}
	nBarcodesEstimate <- replicate(n = 100, estimateBarcodes() )

	results <- rbind(results, data.frame(depth = sample_depth, barcodes = nBarcodesEstimate))
	
}

# str(results)

# Plot output
p <- ggplot(results, aes(depth, barcodes))
p <- p + stat_smooth() + geom_point(size = 1)
ggsave(plot = p, path = outputdir, file = "Downsampling_plot.pdf")

