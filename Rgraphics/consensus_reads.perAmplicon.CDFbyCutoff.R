

samplename <- as.character(commandArgs(trailingOnly = TRUE)[1])
PDFfilename <- paste(c(samplename, "ConsensusCount_plots.byAmplicon.pdf"), collapse="_")

# Create Consensus Count composition plots

UIDdata <- read.delim(file=paste(c(samplename, ".consensusSequences.txt.gz"), collapse=""), header=F, 
	colClasses=c("factor", "character", "numeric", "character", "character"))
	
M <- UIDdata[,c(1,3)]
colnames(M) <- c("Amplicon", "Depth")
# M <- M[!(M$Amplicon == "UNKNOWN"),]   # Remove the UNKNOWN amplicons
# str(M)

# Generate the data frame containing the CDF data for plotting
M.cdf <- data.frame(table(M[M$Depth >= 1,"Amplicon"]), d=1)
for ( i in c(2:100, seq(200, 10000, 100)) ) {
	M.cdf <- rbind(M.cdf, data.frame(table(M[M$Depth >= i,"Amplicon"]), d=i))
	}	
# M.cdf

# Create the plot
library(ggplot2)
g <- ggplot(M.cdf, aes(x=d, y=Freq, colour=Var1))
g <- g + geom_line() + scale_y_log10() + scale_x_log10()
g <- g + labs(x = "Family Depth Cutoff", y = "n Consensus Sequences", title = "Consensus", colour = "Amplicon")
ggsave(plot = g, filename=PDFfilename)



# SVNID: $Id: consensus_reads.perAmplicon.CDFbyCutoff.R 319 2015-06-10 14:06:04Z pkrzyzanowski $
