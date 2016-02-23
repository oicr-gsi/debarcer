

mkdir -p testresults
cd testresults
cp /u/pkrzyzanowski/projects/simsenseq/debarcer/demodata/Sample_9297.R1.fastq.gz .

module load debarcer/0.2.0
runDebarcer.sh Sample_9297.R1.fastq.gz Sample_9297

