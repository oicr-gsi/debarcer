
import gzip
import os
import sys
import argparse
import configparser
from itertools import zip_longest


def parse_prep(prepname, prepfile):
    '''
    (str, file) --> configparser.SectionProxy
    Returns parameters for prepname specified in the config file
    '''

    preps = configparser.ConfigParser()
    preps.read(prepfile)
    return preps[prepname.upper()]


def getread(fastq_file):
    """
    (file) -- > itertools.zip_longest
    Takes a fastq file open for reading (in plain text mode) and returns an
    iterator slicing the fastq into 4-line reads
    """
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)


def extract_umis(reads, umi_locs, umi_lens):
    '''
    (list, list, list) -> list
    Take a list of read sequences, a list with the 1-base location of umis
    in the reads (eg umi_locs = [1], umi is located in 1st read of reads, reads[0]),
    a list of umi lengths corresponding to the umi in each read and return a list
    of umi sequences
    '''
    
    # make a list with all umis
    umis = []
    
	 #Iterate through (umi_loc, umi_len) pairs in list of tuple pairs
    for umi_loc, umi_len in zip(umi_locs, umi_lens):
        # get the read with the umi convert 1-base to 0-base position
        read = reads[int(umi_loc) - 1]
        # slice the read to extract the umi sequence
        umis.append(read[0:int(umi_len)])
    return umis


def verify_spacer(reads, umis, spacer_seq):
    '''
    (list, list, str) -> bool
    Take a list of read sequences, a list of umi sequence and a spacer sequence
    string and return True by default or False if spacer not in read immediately after umi 
    '''
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
    Pre-condition: fasqs have the same number of reads and files are in sync
    """
    
    # get the parameters for prepname from the config
    prep = parse_prep(prepname, prepfile)

    # get the number of input (1-3) and reheadered read files (1-2)
    num_reads, actual_reads  = int(prep['INPUT_READS']), int(prep['OUTPUT_READS'])
    # get the indices of reads with  UMI (1-3)
    umi_locs = [int(x.strip()) for x in prep['UMI_LOCS'].split(',')]
    # get the length of the umis (1-100)
    umi_lens = [int(x.strip()) for x in prep['UMI_LENS'].split(',')]
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

    # set default spacer lengths for read1 and read2 files 
    spacer_len_r1, spacer_len_r2 = 0, 0

    # retrieve spacer length if exists
    if spacer:
        spacer_len_r1 = len(spacer_seq)
        # update spacer length for read2 if umi in read2  
        if len(umi_locs) > 1:
            spacer_len_r2 = len(spacer_seq)

    # get the length of the umi for read1 and read2, set to 0 if only in read1
    umi_len_r1 = umi_lens[0]
    if len(umi_lens) > 1:
        umi_len_r2 = umi_lens[1]
    else:
        umi_len_r2 = 0

    print("Preprocessing reads...")


    # check the number of input fastqs
    if num_reads == 3:
        # check the number of output fastqs
        if actual_reads == 2:
            # this assumes that umis are located in a single input file
            # and that no spacer is provided ie. HALOPLEX and SURESELECT
            # create iterators with slices of 4 read lines from each file
            # loop over reads in each file 
            for read1, read2, read3 in zip(getread(r1), getread(r2), getread(r3)):
                # extract umi sequences from read2
                umis = extract_umis([read1[1], read2[1], read3[1]], umi_locs, umi_lens)
                # edit read names from r1 and r3
                read_name1, rest1 = read1[0].rstrip().split(' ')
                read_name2, rest2 = read3[0].rstrip().split(' ')
                # add umi seq to read1 name and write read1 to output file 1
                r1_writer.write(read_name1 + ":" + umis[0] + " " + rest1 + "\n")
                for i in range(1, len(read1)):
                    r1_writer.write(read1[i])
                # add umi seq to read3 name and write read3 to output file 2
                r2_writer.write(read_name2 + ":" + umis[0] + " " + rest2 + "\n")
                for i in range(1, len(read3)):
                    r2_writer.write(read3[i])
        else:
            raise ValueError("Invalid configuration of reads/actual reads.")
    elif num_reads == 2:
        # check if paired end or single end
        if actual_reads == 2:
            # loop over read in each file
            for read1, read2 in zip(getread(r1), getread(r2)):
                # extract umis from read1 and read2
                umis = extract_umis([read1[1], read2[1]], umi_locs, umi_lens)
                
                
                # skip reads without spacer
                if spacer and not verify_spacer([read1[1], read2[1]], umis, spacer_seq):
                    next


                # edit read names from r1 and r2 
                read_name1, rest1 = read1[0].rstrip().split(' ')
                read_name2, rest2 = read2[0].rstrip().split(' ')
                # add umis as a single string to read1 name and write read 1 to output file 1
                r1_writer.write(read_name1 + ":" + ''.join(umis) + " " + rest1 + "\n")
                # remove umi and spacer from read seq. write read to output file 1
                r1_writer.write(read1[1][umi_len_r1 + spacer_len_r1:])
                r1_writer.write(read1[2])
                r1_writer.write(read1[3][umi_len_r1 + spacer_len_r1:])
                # add umis as a single string to read2 name and write read 1 to output file 2 
                r2_writer.write(read_name2 + ":" + ''.join(umis) + " " + rest2 + "\n")
                # remove umi and spacer from read seq. write read to output file 2
                r2_writer.write(read2[1][umi_len_r2 + spacer_len_r2:])
                r2_writer.write(read2[2])
                r2_writer.write(read2[3][umi_len_r2 + spacer_len_r2:])
        # check if paired end or single end
        elif actual_reads == 1:
            # loop over reads in each file
            for read1, read2 in zip(getread(r1), getread(r2)):
                # extract umis from read1 and read2
                umis = extract_umis([read1[1], read2[1]], umi_locs, umi_lens)
                # edit read names from r1
                read_name1, rest1 = read1[0].rstrip().split(' ')
                # add # add umi seq to read1 name and write read1 to output file 1
                r1_writer.write(read_name1 + ":" + umis[0] + " " + rest1 + "\n")
                # write remaining of read1 to output file 1
                r1_writer.write(read1[1])
                r1_writer.write(read1[2])
                r1_writer.write(read1[3])
        else:
            raise ValueError("Invalid configuration of reads/actual reads.")
    else:
        if actual_reads == 1:
            # loop over reads in r1
            for read1 in getread(r1):
                # extract umi from read1
                umis = extract_umis([read1[1]], umi_locs, umi_lens)
                # skip parser without reads
                if spacer and not verify_spacer([read1[1], read2[1]], umis, spacer_seq):
                    next

                # edit read name
                read_name1, rest1 = read1[0].rstrip().split(' ')
                # add umi to read name and write to outputfile
                r1_writer.write(read_name1 + ":" + umis[0] + " " + rest1 + "\n")
                # remove umi and spacer from read seq. write remaining of read to output file
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

    print("Complete. Output written to {}.".format(outdir))


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
    