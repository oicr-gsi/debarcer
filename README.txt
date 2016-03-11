

Debarcer: De-Barcoding and Error Correction

You can export or checkout the whole package from the following location:

  svn export file:///u/pkrzyzanowski/svn/repo/projects/EAC/molecular_barcoding/trunk/debarcer .
  svn co file:///u/pkrzyzanowski/svn/repo/projects/EAC/molecular_barcoding/trunk/debarcer .
  
  

Running Debarcer  
----------------

To generate a set of positional compositions, you can either export 
debarcer to the location, or only the exec script:

  svn export file:///u/pkrzyzanowski/svn/repo/projects/EAC/molecular_barcoding/trunk/Debarcer/generatePositionalCompositionFiles.sh

--> New (Oct 2014):
You can execute the script directly from the Debarcer directory:

  /u/pkrzyzanowski/projects/EAC/molecular_barcoding/svn/trunk/Debarcer/generatePositionalCompositionFiles.sh Sample_6.R1.fastq.gz Sample6

Or via qsub:

  qsub -N "BCS" -b y -cwd -l h_vmem=16g "/u/pkrzyzanowski/projects/EAC/molecular_barcoding/svn/trunk/Debarcer/runDebarcer.sh Sample_1_ad4.R1.fastq.gz Sample01"

To use alternate tables (5-plex vs 15-plex), you need to copy the script over and edit 
the AMPLICON_TABLE variable.

Requirements:  Consensus calling script needs at least 16G memory.

Miscellaneous
-------------

For testing, export the Boarat directory to the environment

  export DEBARCER='/u/pkrzyzanowski/projects/EAC/molecular_barcoding/svn/trunk/Debarcer'


