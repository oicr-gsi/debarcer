# Run this in the root of a results hierarchy to 
# create a tar.gz containing all the Debarcer output
# files
#
# Paul Krzyzanowski
# pmkrzyzanowski@gmail.com
# (c) 2016


STARTTIME=$(date +"%Y%m%d-%H%M%S")
RESULTSFILENAME="Debarcer_Results_"$STARTTIME".tar.gz"

tar cvfz $RESULTSFILENAME `find . -name *.pdf \
	-o -name "*bamPositionComposition*" \
	-o -name "*UID*" \
	-o -name "*SummaryStatistics.txt" \
	-o -name "*log"`


