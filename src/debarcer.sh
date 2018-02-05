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
    -d|--bed_file)
    BED_FILE="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--python)
    PYTHON="$2"
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

# Load Python from specified script (revisit, TODO)
bash $PYTHON

# 1. Get bed regions
python get_bed_regions.py "$BED_FILE" "$OUTPUT" "$REGION"

# 2. Perform UMI tally
python UMI_count.py "$OUTPUT" "$BAM_FILE" "$OUTPUT/$CHR:$POS_A-$POS_B.regions" "$REGION"

# 3. Generate consensus
python generate_consensus.py "$BAM_FILE" "$OUTPUT/$CHR:$POS_A-$POS_B.tally" "$CONFIG" "$REGION"

# 4. TODO stats/plots/etc from consensus...
python generate_report.py "$OUTPUT/$CHR:$POS_A-$POS_B.cons" "$OUTPUT" "$BAM_FILE" "$REGION" "$CONFIG"

