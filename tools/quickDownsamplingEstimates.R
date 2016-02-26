


UIDfilename <- "Sample_9297.UIDdepths.txt.gz"
outputdir <- "./downsampling/"
library(ggplot2)

UIDs <- read.delim(file = UIDfilename, header = F)
colnames(UIDs) <- c("Amplicon", "Barcode", "Count")

allAmplicons <- as.character(unique(UIDs$Amplicon))

a <- "chr17:7577046"


UIDs.a <- UIDs[UIDs$Amplicon == a,]

results <- data.frame(depth = 0, barcodes = 0)

for ( d in c(100, 1000, seq(2000, 100000, 2000)) ) {
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

str(results)

p <- ggplot(results, aes(depth, barcodes))
p <- p + stat_smooth() + geom_point()
p
ggsave(path = outputdir, file = "test.pdf")

