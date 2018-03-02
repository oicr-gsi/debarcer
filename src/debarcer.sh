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
    -r1|--read1)
    READ1="$2"
    shift # past argument
    shift # past value
    ;;
    -r2|--read2)
    READ2="$2"
    shift # past argument
    shift # past value
    ;;
    -R|-ref|--reference)
    REF="$2"
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

# 1. TODO Reheader (add option for this)

# 2. Align, sort, index FASTQs if BAM file does not exist
if ! [ -f "$BAM_FILE" ]
then
    module load bwa/0.7.12
    module load samtools
    
    bwa mem -M -t 8 $REF $READ1 $READ2 | samtools view -bS -> "$BAM_FILE.unsorted"
    
    samtools sort "$BAM_FILE.unsorted" "${BAM_FILE/.bam/}"
    samtools index "$BAM_FILE"
fi

# Load Python TODO
# module load python-gsi/3.6.4
# module load pysam
bash $PYTHON 

# 3. Perform UMI tally
python UMI_count.py --bam_file $BAM_FILE --bed_file $BED_FILE --region $REGION --output_path $OUTPUT --config $CONFIG

# 4. Generate consensus
python generate_consensus.py --bam_file $BAM_FILE --tally $OUTPUT/$CHR:$POS_A-$POS_B.tally  --output_path $OUTPUT --region $REGION --config $CONFIG

# 5. TODO stats/plots/etc from consensus...
# python3.6 generate_report.py --tally $OUTPUT/$CHR:$POS_A-$POS_B.fsize1.cons --output_path $OUTPUT --bam_file $BAM_FILE --region $REGION --config $CONFIG

