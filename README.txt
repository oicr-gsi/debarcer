

Debarcer: De-Barcoding and Error Correction



Setting up Debarcer
-------------------

1) Download a debarcer release, i.e.:

  wget https://github.com/oicr-gsi/debarcer/archive/v0.2.0.tar.gz

2) Unpack into a convenient directory, which will become the $BHOME
   enviroment variable.

3) Add $BHOME to your PATH. 

   If you're in an environment that uses
   modules, you can install the appropriate module file from the
   $BHOME/utils/modulefiles directory and run

     module load debarcer/[version number]

   to adjust your environment prior to running the pipeline.


Dependencies
------------

- samtools v0.1.19 or higher
- bwa v0.7.12 or higher
- Bio-SamTools-1.41 (supplied)
-- Bio-SamTools required BioPerl
- R v2.11.0 or higher

The pipeline needs to be run on a maching with at least 16GB memory.
It may run with less but this hasn't been tested extensively.

Running Debarcer
----------------

Once $BHOME is set in your path, you can view debarcer usage like so:

  runDebarcer.sh -u
  
To run an analysis, go to a directory with a barcoded fastq file
and run

  runDebarcer.sh -r -f <infile.fastq.gz> -n <SampleName>

You can use the Sample_Test.R1.fastq.gz file in the $BHOME/test_toole
directory to try out the pipeline.

If you'd like to regenerate only the graphics in one of the analysis
directories, run

  runDebarcer.sh -g -f <infile.fastq.gz> -n <SampleName>
  
Acknowledgements
----------------

Levenshtein C code is a modified version of source available here: http://ideone.com/OL38zi
