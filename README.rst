Debarcer
========

A package for De-Barcoding and Error Correction of sequencing data containing molecular barcodes. For information on getting started with Debarcer, see the `Installation Guide`_.

.. _`Installation Guide`: https://github.com/oicr-gsi/debarcer/wiki/1.-Installation-Guide

Note : The repository for the original release of Debarcer (V0.3.1) as described in Nature Protocols can be found under 
https://github.com/oicr-gsi/debarcer/releases/tag/v0.3.1 
https://github.com/oicr-gsi/debarcer/tree/master-old


Configuration
-------------

A sample config file is provided in ``/debarcer/config/sample_config.ini``, and a sample prepfile is provided in ``/debarcer/config/library_prep_types.ini``. Thresholds and paths that you expect to use multiple times should be stored in the config file. The prepfile contains instructions for how different library preps should be handled, for example:

.. code:: ini

	; One library entry
	[LIBRARY_NAME]
	INPUT_READS=3
	OUTPUT_READS=2
	UMI_LOCS=2
	UMI_LENS=10
	SPACER=FALSE

* INPUT_READS:
	Number of unprocessed fastq files (1-3).
* OUTPUT_READS:
	Number of fastq files after reheadering (1-2).
* UMI_LOCS:
	Comma-separated indices of reads containing a UMI (1-3).
* UMI_LENS:
	Comma-separated lengths of UMIs corresponding to UMI_LOCS (1-100).
* SPACER:
	TRUE if a spacer is present, FALSE otherwise (TRUE/FALSE).
* SPACER_SEQ (optional):
	Base sequence of the spacer ([A,C,G,T]+).

``LIBRARY_NAME`` would then be the ``--prepname`` argument passed to debarcer.py.

Typical Workflow
----------------

.. code:: bash

	## Preprocess some fastq files
	$ python debarcer.py preprocess -o /path/to/output_dir -r1 /path/to/read1.fastq -r2 /path/to/read2.fastq
	  -prepname "prepname" -prepfile /path/to/library_prep_types.ini

	## Align, sort, index
	## ...
	## produces: bam_file.bam, bam_file.bam.bai

	## Error-correct and group UMIs into families
	$ python debarcer.py group -r chrN:posA-posB -c /path/to/config.ini -b /path/to/bam_file.bam
	  -o /path/to/output_dir

	## Perform base collapsing
	$ python debarcer.py collapse -o /path/to/output_dir -r chrN:posA-posB
          -b /path/to/bam_file.bam -u /path/to/umi_file.umis
          -c /path/to/config.ini
	  
	## Call variants of specified family sizes
        $ python debarcer.py call -o /path/to/output_dir -r chrN:posA-posB
          -cf /path/to/cons_file.cons -f 1,2,5 -c path/to/config.ini

Dependencies
------------

Debarcer was tested using Python 3.6.4 and depends on the packages pysam and pandas. See ``requirements.txt``.
