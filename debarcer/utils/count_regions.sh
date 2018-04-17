BED_FILE=$1
BAM_FILE=$2
CDIR=$3
ERROR=$4
OUTPUT=$5

BED_SCRIPT="$CDIR/get_regions.py"
BAM_SCRIPT="$CDIR/UMI_count.py"

BED_FILE="$CDIR/bed_regions.csv"

while IFS=, read -r contig start end
do    
    qsub -cwd -b y -l h_vmem=32g -e "$ERROR" -o "$OUTPUT/sge" -N "UMI_count" \
        "cd /u/tbodak/python3.6.4-env/bin;source activate;cd $CDIR;python $BAM_SCRIPT $BAM_FILE $OUTPUT $contig $start $end $BED_FILE";

done < regions.csv
