#!/bin/bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -c|--config)
    CONFIG="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--region)
    REGION="$2"
    shift # past argument
    shift # past value
    ;;
    -b|--bam_file)
    BAM_FILE="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

CHR=${REGION%:*}
POS_A=${REGION#*:}
POS_A=${POS_A%-*}
POS_B=${REGION#*-}

# need to find better way to load python
cd /u/tbodak/python3.6.4-env/bin
source activate
cd -

# UMI Count
python UMI_count.py "$OUTPUT" "$BAM_FILE" "$CHR" "$POS_A" "$POS_B"

# consensus
python generate_consensus.py "$BAM_FILE" "$OUTPUT/output_$CHR-$POS_A-$POS_B.txt" "$CONFIG"

# TODO stats/plots/etc from consensus...
python stats.py "$OUTPUT/cons_$CHR-$POS_A-$POS_B.txt" "$OUTPUT" "$CONFIG"

