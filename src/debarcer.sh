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

# 1. Perform UMI tally
python UMI_count.py --bam_file $BAM_FILE --bed_file $BED_FILE --region $REGION --output_path $OUTPUT --config $CONFIG

# 2. Generate consensus
python generate_consensus.py --bam_file $BAM_FILE --tally $OUTPUT/$CHR:$POS_A-$POS_B.tally  --output_path $OUTPUT --region $REGION --config $CONFIG

# 3. TODO stats/plots/etc from consensus...
# python generate_report.py --tally $OUTPUT/$CHR:$POS_A-$POS_B.fsize1.cons --output_path $OUTPUT --bam_file $BAM_FILE --region $REGION --config $CONFIG

