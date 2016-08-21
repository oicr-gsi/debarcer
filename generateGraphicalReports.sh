

# Generate the graphics from files created by the main runDebarcer.sh script
#

DBROOT=$1
SAMPLENAME=$2
GRAPHICSSRC=$DBROOT"/src/graphics"

# Plotting functions
module load R/3.3.0

# Rscript $GRAPHICSSRC/analysis_plotting.R.txt $SAMPLENAME  # For ecdfs, FIXME

# Generate the error rate dotplots
Rscript $GRAPHICSSRC/bam_analysis_plotting.barplot.txt $SAMPLENAME cons3 # Need to supply the suffix for a bamPositionComposition.???.txt file
Rscript $GRAPHICSSRC/bam_analysis_plotting.barplot.txt $SAMPLENAME cons10 # Need to supply the suffix for a bamPositionComposition.???.txt file
Rscript $GRAPHICSSRC/bam_analysis_plotting.barplot.txt $SAMPLENAME cons20 # Need to supply the suffix for a bamPositionComposition.???.txt file
Rscript $GRAPHICSSRC/bam_analysis_plotting.barplot.txt $SAMPLENAME cons30 # Need to supply the suffix for a bamPositionComposition.???.txt file

# These are the dodged position by position error comparisons bar plots
Rscript $GRAPHICSSRC/bam_analysis.ErrorBarplots.txt $SAMPLENAME cons3 # Need to supply the suffix for a bamPositionComposition.???.txt file
Rscript $GRAPHICSSRC/bam_analysis.ErrorBarplots.txt $SAMPLENAME cons10 # Need to supply the suffix for a bamPositionComposition.???.txt file
Rscript $GRAPHICSSRC/bam_analysis.ErrorBarplots.txt $SAMPLENAME cons20 # Need to supply the suffix for a bamPositionComposition.???.txt file
Rscript $GRAPHICSSRC/bam_analysis.ErrorBarplots.txt $SAMPLENAME cons30 # Need to supply the suffix for a bamPositionComposition.???.txt file

# Rscript $GRAPHICSSRC/analysis_plotting.barplot.txt $SAMPLENAME bc20 # Need to supply the suffix for a consensusPositionalComposition.???.txt file
# Rscript $GRAPHICSSRC/analysis_plotting.barplot.txt $SAMPLENAME bc30 # Need to supply the suffix for a consensusPositionalComposition.???.txt file
# Rscript $GRAPHICSSRC/analysis_plotting.barplot.txt $SAMPLENAME bc100 # Need to supply the suffix for a consensusPositionalComposition.???.txt file
# Rscript $GRAPHICSSRC/analysis_plotting.AmpliconsByPosition.txt $SAMPLENAME bc3  
# Rscript $GRAPHICSSRC/analysis_plotting.AmpliconsByPosition.txt $SAMPLENAME bc10
# Rscript $GRAPHICSSRC/analysis_plotting.AmpliconsByPosition.txt $SAMPLENAME bc20
# Rscript $GRAPHICSSRC/analysis_plotting.AmpliconsByPosition.txt $SAMPLENAME bc30
# Rscript $GRAPHICSSRC/analysis_plotting.R.bc10strict.txt $SAMPLENAME

Rscript $GRAPHICSSRC/UID_depth_plotting.R.txt $SAMPLENAME
Rscript $GRAPHICSSRC/UID_depth_plotting.perAmplicon.R.txt $SAMPLENAME  # Generate PDF with many plots of UID/depth distributions, per amplicon
# Rscript $GRAPHICSSRC/consensus_reads.perAmplicon.CDFbyCutoff.R $SAMPLENAME  # Generate curves of consensus read counts across a range of minimum depths
