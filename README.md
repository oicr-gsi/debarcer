# debarcer
Debarcer: A package for De-Barcoding and Error Correction of sequencing data containing molecular barcodes

Typical Workflow
-----
```shell
## Preprocess some fastq files
python debarcer.py P -o /path/to/output_dir -r1 /path/to/read1.fastq -r2 /path/to/read2.fastq -prepname "prepname" -prepfile /path/to/library_prep_types.ini

## Align, sort, index 
## ...
## produces: bam_file.bam, bam_file.bam.bai

## Run Debarcer 
python debarcer.py D -t -r chrN:posA-posB -c /path/to/config.ini  -be /path/to/bed_file.bed -b /path/to/bam_file.bam -o /path/to/output_dir
```

A sample config file is provided in `/debarcer/config/sample_config.ini`, and a sample prepfile is provided in `/debarcer/config/library_prep_types.ini`.
