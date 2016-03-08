

mkdir -p testresults
cd testresults

module load debarcer/dev
cp $BHOME/test_tools/Sample_Test.R1.fastq.gz .
runDebarcer.sh Sample_Test.R1.fastq.gz Sample_Test