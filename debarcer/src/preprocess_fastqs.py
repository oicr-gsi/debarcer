import gzip
import os
import argparse
import configparser
from itertools import zip_longest


def parse_prep(prepname, prepfile):
    '''
    (str, file) --> configparser.SectionProxy
    :param prepname: Name of the library preparation
    :param prepfile: Path to the library preparation ini file
    
    Returns key, value pairs of parameters for prepname specified in the prepfile
    '''

    preps = configparser.ConfigParser()
    preps.read(prepfile)
    return preps[prepname.upper()]


def getread(fastq_file):
    """
    (file) -- > itertools.zip_longest
    :param fastq_file: a fastq file open for reading in plain text mode
    
    Returns an iterator slicing the fastq into 4-line reads.
    Each element of the iterator is a tuple containing read information
    """
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)


def extract_umis(reads, umi_locs, umi_lens):
    '''
    (list, list, list) -> list
    :param reads: A list of read sequences
    :param umi_locs: A list of 1-based indices indicating which read sequences have the umis.
     (eg umi_locs = [1]: umi is located in 1st read of reads, reads[0])
    :param umi_lens: A list with umi lengths for each location
    
    Returns a return a list of umi sequences
    '''
    
    # make a list with all umi sequences
    umis = []
    
	 #Iterate through (umi_loc, umi_len) pairs in list of tuple pairs
    for umi_loc, umi_len in zip(umi_locs, umi_lens):
        # get the read with the umi convert 1-base to 0-base position
        read = reads[int(umi_loc) - 1]
        # slice the read to extract the umi sequence
        umis.append(read[0:int(umi_len)])
    return umis


def correct_spacer(reads, umis, spacer_seq):
    '''
    (list, list, str) -> bool
    Take a list of read sequences, a list of umi sequence and a spacer sequence
    string and return True by default or False if spacer not in read immediately after umi 
    '''
    
    # set up bool <- True
    # by default: spacer seq follows umi seq in read sequences with spacer 
    Correct = True
    # for each read seq and umi seq
    #     Correct <- False if umi in read and spacer position not as expected 
    for read in reads:
        for umi in umis:
            # check if umi in read. ignore reads if umi not in read
            if umi in read:
                # remove umi from read. 
                read = read[read.index(umi) + len(umi):]
                # check is spacer seq is immediately after umi 
                if not read.upper().startswith(spacer_seq.upper()):
                    # update bool
                    Correct = False
    return Correct


def reheader_fastqs(r1_file, r2_file, r3_file, outdir, prefix, prepname, prepfile):
    """
    (str, str, str, str, str, str, str) -> None
    Take at least 1 input fastq file (r1_file), and reheader fastq file(s)
    according to the prename-specified library prep. in prepfile. 
    Reheadered fastqs are written in outdir and named prefix.umi.reheadered_RN.fastq.gz       
    - removes reads without a valid spacer (if applicable)
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
    
    # get the spacer sequence and spacer length if exists 
    if spacer:
        spacer_seq = str(prep['SPACER_SEQ'])
        spacer_len_r1 = len(spacer_seq)
        # update spacer length for read2 if umi in read2  
        if len(umi_locs) > 1:
            spacer_len_r2 = len(spacer_seq)
        else:
            spacer_len_r2 = 0
    else:
        spacer_seq = 'None'
        spacer_len_r1, spacer_len_r2 = 0, 0    
    
    # Read FASTQ in text mode
    r1 = gzip.open(r1_file, "rt")
    r2 = gzip.open(r2_file, "rt") if num_reads > 1 else None
    r3 = gzip.open(r3_file, "rt") if num_reads > 2 else None
    # Open output fastqs in text mode for writing
    r1_writer = gzip.open(os.path.join(outdir, prefix + ".umi.reheadered_R1.fastq.gz"), "wt")
    r2_writer = gzip.open(os.path.join(outdir, prefix + ".umi.reheadered_R2.fastq.gz"), "wt") if actual_reads > 1 else None

    # get the length of the umi for read1 and read2, set to 0 if only in read1
    umi_len_r1 = umi_lens[0]
    if len(umi_lens) > 1:
        umi_len_r2 = umi_lens[1]
    else:
        umi_len_r2 = 0

    # do a check based on number of input and output files currently supported in the library prep ini
    if num_reads == 3:
        assert actual_reads == 2, 'Expecting 2 output fastqs and 3 input fastqs'
    elif num_reads == 2:
        assert actual_reads == 2, 'Expecting 2 output fastqs and 2 input fastqs'
    elif num_reads == 1:
        assert actual_reads == 1, 'Expecting 1 output fastq and 1 input fastq'

    print("Preprocessing reads...")
    
#    # check the number of input fastqs
#    if num_reads == 3:
#        # check the number of output fastqs
#        if actual_reads == 2:
#            # loop over iterator with slices of 4 read lines from each file
#            for read1, read2, read3 in zip(getread(r1), getread(r2), getread(r3)):
#                # extract umi sequences from read2
#                umis = extract_umis([read1[1], read2[1], read3[1]], umi_locs, umi_lens)
#                
#                # skip reads with spacer in wrong position
#                if spacer == True and correct_spacer([read1[1], read2[1], read3[1]], umis, spacer_seq) == False:
#                    continue
#                
#                # edit read names from r1 and r3
#                read_name1, rest1 = read1[0].rstrip().split(' ')
#                read_name2, rest2 = read3[0].rstrip().split(' ')
#                # add umi seq to read1 name and write read1 to output file 1
#                r1_writer.write(read_name1 + ":" + ''.join(umis) + " " + rest1 + "\n")
#                # remove umi and spacer from if present and write read to output file 1
#                r1_writer.write(read1[1][umi_len_r1 + spacer_len_r1:])
#                r1_writer.write(read1[2])
#                r1_writer.write(read1[3][umi_len_r1 + spacer_len_r1:])
#                # add umi seq to read3 name and write read3 to output file 2
#                r2_writer.write(read_name2 + ":" + ''.join(umis) + " " + rest2 + "\n")
#                # remove umi and spacer from if present and write read to output file 2
#                r2_writer.write(read3[1][umi_len_r2 + spacer_len_r2:])
#                r2_writer.write(read3[2])
#                r2_writer.write(read3[3][umi_len_r2 + spacer_len_r2:])
#        else:
#            raise ValueError("Invalid configuration of reads/actual reads.")
#    elif num_reads == 2:
#        if actual_reads == 2:
#            # configuration assumes spacer. paired end reads. (ie EPIC-DS and SIMSENSEQ-PE)
#            # oop over iterator with slices of 4 read lines from each file
#            for read1, read2 in zip(getread(r1), getread(r2)):
#                # extract umis from read1 and read2
#                umis = extract_umis([read1[1], read2[1]], umi_locs, umi_lens)
#                
#                # skip reads with spacer in wrong position
#                if spacer == True and correct_spacer([read1[1], read2[1]], umis, spacer_seq) == False:
#                    continue
#                
#                # edit read names from r1 and r2 
#                read_name1, rest1 = read1[0].rstrip().split(' ')
#                read_name2, rest2 = read2[0].rstrip().split(' ')
#                # add umis as a single string to read1 name and write read 1 to output file 1
#                r1_writer.write(read_name1 + ":" + ''.join(umis) + " " + rest1 + "\n")
#                # remove umi and spacer from read seq. write read to output file 1
#                r1_writer.write(read1[1][umi_len_r1 + spacer_len_r1:])
#                r1_writer.write(read1[2])
#                r1_writer.write(read1[3][umi_len_r1 + spacer_len_r1:])
#                # add umis as a single string to read2 name and write read 1 to output file 2 
#                r2_writer.write(read_name2 + ":" + ''.join(umis) + " " + rest2 + "\n")
#                # remove umi and spacer from read seq if umi and spacer in r2. write read to output file 2
#                r2_writer.write(read2[1][umi_len_r2 + spacer_len_r2:])
#                r2_writer.write(read2[2])
#                r2_writer.write(read2[3][umi_len_r2 + spacer_len_r2:])
#        else:
#            raise ValueError("Invalid configuration of reads/actual reads.")
#    elif num_reads == 1:
#        if actual_reads == 1:
#            # loop over iterator with slices of 4 read lines from file r1
#            for read1 in getread(r1):
#                # extract umi from read1
#                umis = extract_umis([read1[1]], umi_locs, umi_lens)
#                
#                # skip reads with spacer in wrong position
#                if spacer == True and correct_spacer([read1[1]], umis, spacer_seq) == False:
#                    continue
#
#                # edit read name
#                read_name1, rest1 = read1[0].rstrip().split(' ')
#                # add umi to read name and write to outputfile
#                r1_writer.write(read_name1 + ":" + umis[0] + " " + rest1 + "\n")
#                # remove umi and spacer from read seq. write remaining of read to output file
#                r1_writer.write(read1[1][umi_len_r1 + spacer_len_r1:])
#                r1_writer.write(read1[2])
#                r1_writer.write(read1[3][umi_len_r1 + spacer_len_r1:])
#        else:
#            raise ValueError("Invalid configuration of reads/actual reads.")

    # make a list of fastqs open for reading
    fastqs = [i for i in [r1, r2, r3] if i != None]
    
    # make a list of files open for writing
    writers = [i for i in [r1_writer, r2_writer] if i != None]
    
    # check the number of input files
    # create iterator with reads from each file
    if len(fastqs) == 3:
        I =  zip(getread(fastqs[0]), getread(fastqs[1]), getread(fastqs[2]))
    elif len(fastqs) == 2:
        I =  zip(getread(fastqs[0]), getread(fastqs[1]))
    elif len(fastqs) == 1:
        I =  zip(getread(fastqs[0]))
        
    # loop over iterator with slices of 4 read lines from each line
    for reads in I:
        # extract umi sequences from reads
        # make a list of read sequences
        readseqs = [i[1] for i in reads]
        umis = extract_umis(readseqs, umi_locs, umi_lens)
        # skip reads with spacer in wrong position
        if spacer == True and correct_spacer(readseqs, umis, spacer_seq) == False:
            continue
        
        # edit read names and add umi
        # make parallel lists with begining and end of read name from r1, and from r2 or r3
        read_name1, rest1 = reads[0][0].rstrip().split(' ')
        readnames, namerests = [read_name1], [rest1]
        
        if len(reads) > 1:
            # edit read name from r2 (or r3)
            read_name2, rest2 = reads[-1][0].rstrip().split(' ')
            readnames.append(read_name2)
            namerests.append(rest2)
         
        # make lists with umi lengths and spacer lengths    
        UmiLength, SpacerLength = [umi_len_r1, umi_len_r2], [spacer_len_r1, spacer_len_r2]    
        
        for i in range(len(writers)):
            # add umi to read name and write to outputfile
            writers[i].write(readnames[i] + ":" + ''.join(umis) + " " + namerests[i] + "\n")
            # remove umi and spacer from read seq. write reamining of read to outputfile
            if i == 0:
                # determine index for reads k <- 0 for r1, -1 for r2 or r3
                k = 0
            elif i > 0:
                k = -1
            writers[i].write(reads[k][1][UmiLength[i] + SpacerLength[i]:])
            writers[i].write(reads[k][2])
            writers[i].write(reads[k][3][UmiLength[i] + SpacerLength[i]:])
            
    # close all open files
    for i in writers:
        i.close()
    for i in fastqs:
        i.close()
        
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
    