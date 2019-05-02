
import gzip
import os
import sys
import argparse
import configparser
from itertools import zip_longest


def parse_prep(prepname, prepfile):
    '''
    (str, file) --> configparser.SectionProxy
    Return parameters for prepname specified in the config file
    '''

    preps = configparser.ConfigParser()
    preps.read(prepfile)
    return preps[prepname.upper()]


def getread(fastq_file):
    """(Iter) slices fastq into 4-line reads."""
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)


def extract_umis(reads, umi_locs, umi_lens):
    """Gets the UMI from a read."""
    umis = []

	#Iterate through (umi_loc, umi_len) pairs in list of tuple pairs
    for umi_loc, umi_len in zip(umi_locs, umi_lens):
        read = reads[int(umi_loc) - 1]
        umis.append(read[0:int(umi_len)])

    return umis


def verify_spacer(reads, umis, spacer_seq):
    """Returns true if spacer is present."""
    for read in reads:
        for umi in umis:
            if umi in read:
                if not read.split(umi)[1].startswith(spacer_seq):
                    return False
    return True


def reheader_fastqs(r1_file, r2_file, r3_file, outdir, prefix, prepname, prepfile):
    """
    (str, str, str, str, str, str, str) -> None
    Take at least 1 input fastq file (r1_file), and reheader fastq file(s)
    according to the prename-specified library prep. in prepfile. 
    Reheadered fastqs are written in outdir and named prefix.umi.reheadered_RN.fastq.gz       
    - removes reads without a valid spacer (if applicable)
    - gzip module is very slow, consider subprocess (at the cost of compatibility)
    """
    
    # get the parameters for prepname from the config
    prep = parse_prep(prepname, prepfile)

    # get the number of input (1-3) and reheadered read files (1-2)
    num_reads, actual_reads  = int(prep['INPUT_READS']), int(prep['OUTPUT_READS'])
    # get the indices of reads with  UMI (1-3)
    umi_locs = [str(x) for x in prep['UMI_LOCS'].split(',')]
    # get the length of the umis (1-100)
    umi_lens = [int(x) for x in prep['UMI_LENS'].split(',')]
    # specify if a spacer is used or not
    spacer = bool(prep['SPACER'])
    
    # get the spacer sequence if exists 
    if spacer:
        spacer_seq = str(prep['SPACER_SEQ'])
    else:
        spacer_seq = None
    
	 # Read FASTQ in text mode
    r1 = gzip.open(r1_file, "rt")
    r2 = gzip.open(r2_file, "rt") if num_reads > 1 else None
    r3 = gzip.open(r3_file, "rt") if num_reads > 2 else None
    # Open output fastqs in text mode for writing
    r1_writer = gzip.open(os.path.join(outdir, prefix + ".umi.reheadered_R1.fastq.gz"), "wt")
    r2_writer = gzip.open(os.path.join(outdir, prefix + ".umi.reheadered_R2.fastq.gz"), "wt") if actual_reads > 1 else None


    spacer_len_r1 = 0
    spacer_len_r2 = 0

    if spacer:
        spacer_len_r1 = len(spacer_seq)

        if len(umi_locs) > 1:
            spacer_len_r2 = len(spacer_seq)

    umi_len_r1 = umi_lens[0]

    if len(umi_lens) > 1:
        umi_len_r2 = umi_lens[1]
    else:
        umi_len_r2 = 0

    print("Preprocessing reads...")

    if num_reads == 3:

        if actual_reads == 2:

            for read1, read2, read3 in zip(getread(r1), getread(r2), getread(r3)):

                umis = extract_umis(
                    [read1[1], read2[1], read3[1]], umi_locs, umi_lens)

                read_name1, rest1 = read1[0].rstrip().split(' ')
                read_name2, rest2 = read3[0].rstrip().split(' ')

                r1_writer.write(read_name1 + ":" +
                                umis[0] + " " + rest1 + "\n")
                r1_writer.write(read1[1])
                r1_writer.write(read1[2])
                r1_writer.write(read1[3])

                r2_writer.write(read_name2 + ":" +
                                umis[0] + " " + rest2 + "\n")
                r2_writer.write(read3[1])
                r2_writer.write(read3[2])
                r2_writer.write(read3[3])

        else:
            raise ValueError("Invalid configuration of reads/actual reads.")

    elif num_reads == 2:

        if actual_reads == 2:

            for read1, read2 in zip(getread(r1), getread(r2)):

                umis = extract_umis([read1[1], read2[1]], umi_locs, umi_lens)

                if spacer and not verify_spacer([read1[1], read2[1]], umis, spacer_seq):
                    next

                read_name1, rest1 = read1[0].rstrip().split(' ')
                read_name2, rest2 = read2[0].rstrip().split(' ')

                r1_writer.write(read_name1 + ":" +
                                ''.join(umis) + " " + rest1 + "\n")
                r1_writer.write(read1[1][umi_len_r1 + spacer_len_r1:])
                r1_writer.write(read1[2])
                r1_writer.write(read1[3][umi_len_r1 + spacer_len_r1:])

                r2_writer.write(read_name2 + ":" +
                                ''.join(umis) + " " + rest2 + "\n")
                r2_writer.write(read2[1][umi_len_r2 + spacer_len_r2:])
                r2_writer.write(read2[2])
                r2_writer.write(read2[3][umi_len_r2 + spacer_len_r2:])

        elif actual_reads == 1:

            for read1, read2 in zip(getread(r1), getread(r2)):

                umis = extract_umis([read1[1], read2[1]], umi_locs, umi_lens)
                read_name1, rest1 = read1[0].rstrip().split(' ')
                r1_writer.write(read_name1 + ":" +
                                umis[0] + " " + rest1 + "\n")
                r1_writer.write(read1[1])
                r1_writer.write(read1[2])
                r1_writer.write(read1[3])

        else:
            raise ValueError("Invalid configuration of reads/actual reads.")

    else:

        if actual_reads == 1:

            for read1 in getread(r1):

                umis = extract_umis([read1[1]], umi_locs, umi_lens)

                if spacer:

                    if not verify_spacer([read1[1], read2[1]], umis, spacer_seq):
                        next

                    spacer_len = len(spacer_seq)

                else:
                    spacer_len = 0

                read_name1, rest1 = read1[0].rstrip().split(' ')
                r1_writer.write(read_name1 + ":" +
                                umis[0] + " " + rest1 + "\n")
                r1_writer.write(read1[1][umi_len_r1 + spacer_len_r1:])
                r1_writer.write(read1[2])
                r1_writer.write(read1[3][umi_len_r1 + spacer_len_r1:])

        else:
            raise ValueError("Invalid configuration of reads/actual reads.")

    r1.close()
    if r2:
        r2.close() 
    if r3:
        r3.close() 

    r1_writer.close()
    if r2_writer:
        r2_writer.close() 

    print("Complete. Output written to {}.".format(output_path))


if __name__ == '__main__':

    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1', '--Read1', dest='read1', help='Path to first FASTQ file.', required=True)
    parser.add_argument('-r2', '--Read2', dest='read2', help='Path to second FASTQ file, if applicable')
    parser.add_argument('-r3', '--Read3', dest='read3', help='Path to third FASTQ file, if applicable')
    parser.add_argument('-p',  '--Prepname', dest='prepname', choices=['HALOPLEX', 'SURESELECT', 'EPIC-DS', 'SIMSENSEQ-PE', 'SIMSENSEQ-SE'],
                        help='Name of library prep to  use (defined in library_prep_types.ini)', required=True)
    parser.add_argument('-pf', '--Prepfile', dest='prepfile', help='Path to the library_prep_types.ini file', required=True)
    parser.add_argument('-o', '--OutDir', dest='outdir', help='Output directory where fastqs are written', required=True)
    parser.add_argument('-px', '--Prefix', dest= 'prefix', help='Prefix for naming umi-reheradered fastqs. Use Prefix from Read1 if not provided') 
    
    args = parser.parse_args()

    r1_file, r2_file, r3_file = args.read1, args.read2, args.read3
    outdir, prefix = args.outdir, args.prefix
    prepname, prepfile = args.prepname, args.prepfile

    # Preprocess (reheader fastq files)
    reheader_fastqs(r1_file, r2_file, r3_file, outdir, prefix, prepname, prepfile)
    