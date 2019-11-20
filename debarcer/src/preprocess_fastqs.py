import gzip
import os
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


def check_library_prep(prepname, prepfile):
    '''
    (str, file) --> None
    :param prepname: Name of the library preparation
    :param prepfile: Path to the library preparation ini file
    
    Check that library prep ini has valid values. Raise ValueError if not
    '''
    
    # get the parameters from the prep file for a given library prep
    L = parse_prep(prepname, prepfile)
    D = {i:j for i, j in L.items()} 
    # check that all expected parameters are present
    expected = {'input_reads', 'output_reads', 'umi_locs', 'umi_lens', 'spacer', 'spacer_seq', 'umi_pos'}
    
    unexpected = set(D.keys()).difference(set(expected))
    if len(unexpected) != 0:
        raise ValueError('ERR: unexpected keys in librabry prep')
    missing = set(expected).symmetric_difference(set(D.keys()))
    if len(missing) != 0:
        raise ValueError('ERR: missing keys in librabry prep')
        
    # check that umi_locs, umi_lens and umi_pos have same length
    a = list(map(lambda x: x.strip(), D['umi_locs'].split(',')))
    b = list(map(lambda x: x.strip(), D['umi_lens'].split(',')))
    c = list(map(lambda x: x.strip(), D['umi_pos'].split(',')))
    
    # if , not in umi_pos or umi_lens, the same values is propagated to all files having umis
    if ',' not in D['umi_lens']:
        b = b * len(a)
    if ',' not in D['umi_pos']:
        c = c * len(a)
    if len(a) != len(b) != len(c):
        raise ValueError('ERR: umi_locs and umi_lens should be comma-separated lists of identical size')
       
    # accepted nucleotides in spacer sequences
    nucleotides = 'ACGTURYSWKMBDHVN'
    # check parameter format    
    for i in D:
        if i in ['input_reads', 'output_reads']:
            try:
                int(D[i])
            except:
                raise ValueError('ERR: value for {0} should be an integer'.format(i))
        elif i in ['umi_locs', 'umi_lens', 'umi_pos']:
            for j in list(map(lambda x: x.strip(), D[i].split(','))):
                try:
                    int(j)
                except:
                    raise ValueError('ERR: value for {0} should be a comma separated list of integers'.format(i))
        elif i == 'spacer':
            if D['spacer'].upper() == 'TRUE':
                D['spacer'] = True
            elif D['spacer'].upper() == 'FALSE':
                D['spacer'] = False
            else:
                raise ValueError('ERR: value for {0} should be boolean'.format(i))
            if D['spacer'] == True:
                if D['spacer_seq'].lower() in ['none', '']:
                    raise ValueError('ERR: spacer_seq if not defined')
                non_valid = set(D['spacer_seq'].upper()).difference(set(nucleotides.upper()))    
                if len(non_valid) != 0:
                    raise ValueError('ERR: spacer sequence contains non valid nucleotides: {0}'.format(', '.join(non_valid)))
                
def getread(fastq_file):
    """
    (file) -- > itertools.zip_longest
    :param fastq_file: a fastq file open for reading in plain text mode
    
    Returns an iterator slicing the fastq into 4-line reads.
    Each element of the iterator is a tuple containing read information
    """
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)


def extract_umis(reads, umi_locs, umi_lens, umi_pos):
    '''
    (list, list, list, list) -> list
    :param reads: A list of read sequences
    :param umi_locs: A list of 1-based indices indicating which read sequences have the umis.
     (eg umi_locs = [1]: umi is located in 1st read of reads, reads[0])
    :param umi_lens: A list with umi lengths for each location
    :param umi_pos: A list of umi_positions (1-based) within reads
    
    
    Returns a return a list of umi sequences
    Pre-condition: umi_locs, umi_lens and umi_pos have same length
    '''
    
    # make a list with all umi sequences
    umis = []
    
	 #Iterate through (umi_loc, umi_len) pairs in list of tuple pairs
    for umi_loc, umi_len, pos in zip(umi_locs, umi_lens, umi_pos):
        # get the read with the umi convert 1-base to 0-base position
        read = reads[int(umi_loc) - 1]
        # convert the left-most umi position to 0-based
        pos = pos - 1
        # slice the read to extract the umi sequence
        umis.append(read[pos: pos + int(umi_len)])
    return umis


def correct_spacer(read, umi, umi_pos, spacer_seq):
    '''
    (str, str, int, str) -> bool
    :param reads: Read sequence
    :param umis: Umi sequence
    :param umi_pos: Expected position 0-based on umi in read
    :spacer_seq: Spacer sequence
    
    Return False if umi not in read or spacer in sequence but not at expected position and return True otherwise 
    '''
    
    # set up bool <- True
    # by default: umi/spacer config is correct 
    Correct = True
    if umi.upper() not in read.upper():
        Correct = False
    else:
        read = read[umi_pos + len(umi):]
        # check is spacer seq is immediately after umi 
        if not read.upper().startswith(spacer_seq.upper()):
            # update bool
            Correct = False
    return Correct

def open_optional_file(D, k):
    '''
    (dict, str) -> None or opened file in plain text reading mode 
    :param D: a dictionary with key, value pairs
    :param k: a specific key, may or may not be in D
    
    Return a file opened for reading in plain text if k in D and D[k] == True or None otherwise
    '''
    
    # check if key k in D
    if k in D:
        if D[k]:
            r = gzip.open(D[k], "rt")
        else:
            r = None
    else:
        r = None
    return r


def extract_prefix_from_filename(fastqfile):
    '''
    (file) -> str
    
    :param fastqfile: Path to FASTQ file
    
    Return a prefix extracted from the name of the fastq file
    '''
    
    filename = os.path.basename(fastqfile)
    # remove extension if possible or keep entire file name as prefix
    if 'fastq.gz' in filename:
        prefix = filename[:filename.index('.fastq.gz')]
    elif 'fq.gz' in filename:
        prefix = filename[:filename.index('.fq.gz')]
    else:
        prefix = filename
    return prefix
    

def reheader_fastqs(r1_file, outdir, prepname, prepfile, **KeyWords):
    """
    (str, str, str, str, str, str, str) -> tuple
    :param r1_file: Path to first FASTQ file
    :param prepname: Name of the library preparation
    :param prepfile: Path to the library preparation ini file
    :param kw: Recognized keyword arguments: \
               prefix: Prefix for naming umi-reheradered fastqs. Use Prefix from r1 if not provided
               r2: Path to second FASTQ file
               r3: Path to third FASTQ file
    
    Write new reheadered fastq file(s) prefix.umi.reheadered_RN.fastq.gz in outdir
    according to settings in prepfile corresponding to prename 
    Returns a tuple with a dictionary with read counts and read length
    (correct, incorrect umi/spacer configuration, total, length read1/read2)
    and a list of all umi sequences (with correct umi/spacer configuration)
    
    Pre-condition: fastqs have the same number of reads and files are in sync
    """
    
    # check if prefix provided   
    if 'prefix' in KeyWords:
        prefix = KeyWords['prefix']
        # extract prefix from r1_file if no value
        if prefix in [None, '']:
            prefix = extract_prefix_from_filename(r1_file)
    else:
        # extract prefix from r1_file
        prefix = extract_prefix_from_filename(r1_file)
        
    # get the parameters for prepname from the config
    prep = parse_prep(prepname, prepfile)

    # get the number of input (1-3) and reheadered read files (1-2)
    num_reads, actual_reads  = int(prep['INPUT_READS']), int(prep['OUTPUT_READS'])
    # get the indices of reads with  UMI (1,2)
    umi_locs = [int(x.strip()) for x in prep['UMI_LOCS'].split(',')]
    # get the length of the umis
    umi_lens = [int(x.strip()) for x in prep['UMI_LENS'].split(',')]
    # get the positions of umis within reads
    umi_pos = [int(x.strip()) for x in prep['UMI_POS'].split(',')]
    #  if a single value is listed in the library_prep.ini for umi_lens and umi_pos, it will be propagated to all fastqs having umis
    if ',' not in prep['UMI_POS']:
        umi_pos = umi_pos *  len(umi_locs)
    if ',' not in prep['UMI_LENS']:
        umi_lens = umi_lens * len(umi_locs)
    
    # specify if a spacer is used or not
    spacer = prep.getboolean('SPACER')
            
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
    
    # open files r2_file and r2_file if provided, return None otherwise
    r2 = open_optional_file(KeyWords, 'r2')
    r3 = open_optional_file(KeyWords, 'r3')
      
    # Open output fastqs in text mode for writing re-headered reads
    r1_writer = gzip.open(os.path.join(outdir, prefix + ".umi.reheadered_R1.fastq.gz"), "wt")
    if r2 != None:
        r2_writer = gzip.open(os.path.join(outdir, prefix + ".umi.reheadered_R2.fastq.gz"), "wt")
    else:
        r2_writer = None

    # write reads with incorrect configuration to new fastqs
    if spacer == True:
        # Open fastqs in text mode for writing discarded reads
        r1_discarded = gzip.open(os.path.join(outdir, prefix + ".incorrect_reads.R1.fastq.gz"), "wt")
        if r2 != None:
            r2_discarded = gzip.open(os.path.join(outdir, prefix + ".incorrect_reads.R2.fastq.gz"), "wt")
        else:
            r2_discarded = None
        if r3 != None:
            r3_discarded = gzip.open(os.path.join(outdir, prefix + ".incorrect_reads.R3.fastq.gz"), "wt")
        else:
            r3_discarded = None
        
    # get the length of the umi for read1 and read2, set to 0 if only in read1
    umi_len_r1 = umi_lens[0]
    if len(umi_lens) > 1:
        umi_len_r2 = umi_lens[1]
    else:
        umi_len_r2 = 0

    # get the umi position 0-based for read1 and read2, set to 0 if only in read1
    # required for removing umi and spacer from read
    # if pos read 2 = 0, spacer length and umi length are also set to 0 -> do not cut read 
    umi_pos_r1 = umi_pos[0] - 1
    if len(umi_pos) > 1:
        umi_pos_r2 = umi_pos[1] - 1
    else:
        umi_pos_r2 = 0

    # do a check based on number of input and output files currently supported in the library prep ini
    if num_reads == 3:
        assert actual_reads == 2, 'Expecting 2 output fastqs and 3 input fastqs'
    elif num_reads == 2:
        assert actual_reads == 2, 'Expecting 2 output fastqs and 2 input fastqs'
    elif num_reads == 1:
        assert actual_reads == 1, 'Expecting 1 output fastq and 1 input fastq'

    print("Preprocessing reads...")
    
    # make a list of fastqs open for reading
    fastqs = [i for i in [r1, r2, r3] if i != None]
    
    # make a list of files open for writing
    writers = [i for i in [r1_writer, r2_writer] if i != None]

    # make a list of files open for writing discarded reads
    if spacer == True:
        discarded = [i for i in [r1_discarded, r2_discarded, r3_discarded] if i != None]
    
    # check the number of input files
    # create iterator with reads from each file
    if len(fastqs) == 3:
        I =  zip(getread(fastqs[0]), getread(fastqs[1]), getread(fastqs[2]))
    elif len(fastqs) == 2:
        I =  zip(getread(fastqs[0]), getread(fastqs[1]))
    elif len(fastqs) == 1:
        I =  zip(getread(fastqs[0]))
    
    # count all reads and reads with incorrect and correct umi/spacer configuration
    Total, Correct, Incorrect = 0, 0, 0
    # Record all umi sequences with correct umi/spacer configuration    
    UmiSequences = []
    # Record read length
    ReadLength = [set(), set()]
        
    # loop over iterator with slices of 4 read lines from each line
    for reads in I:
        # count total reads
        Total += 1
        # extract umi sequences from reads
        # make a list of read sequences
        readseqs = [i[1] for i in reads]
        umis = extract_umis(readseqs, umi_locs, umi_lens, umi_pos)
        
        # make a list of reads with umis
        reads_with_umis = [readseqs[i-1] for i in umi_locs]
        
        # skip reads with spacer in wrong position
        if spacer == True and False in [correct_spacer(reads_with_umis[i], umis[i], umi_pos[i] -1, spacer_seq) for i in range(len(reads_with_umis))]:
            # count reads with incorrect umi/spacer configuration 
            Incorrect += 1
            # write reads with incorrect umi/spacer configuration to separate fastqs 
            for i in range(len(discarded)):
                for j in range(len(reads[i])):
                    discarded[i].write(reads[i][j])
            # skip read processing
            continue

        # count number of reads with correct umi/spacer configuration
        Correct += 1
        
        # collect umi sequences
        UmiSequences.extend(umis)
        
        # edit read names and add umi
        # make parallel lists with begining and end of read name from r1, and from r2 or r3
        read_name1, rest1 = reads[0][0].rstrip().split(' ')
        readnames, namerests = [read_name1], [rest1]
        
        if len(reads) > 1:
            # edit read name from r2 (or r3)
            read_name2, rest2 = reads[-1][0].rstrip().split(' ')
            readnames.append(read_name2)
            namerests.append(rest2)
         
        # make lists with umi lengths, spacer lengths and umi positions for     
        UmiLength, SpacerLength, UmiPositions = [umi_len_r1, umi_len_r2], [spacer_len_r1, spacer_len_r2], [umi_pos_r1, umi_pos_r2]    
        
        for i in range(len(writers)):
            # if paired reads and umis are in each read: concatenate umis and assign concatenated umi to each read
            # if paired read and single umi: assign the same umi to each read
            # if single end read, assign the umi to its read
            if len(umis) > 1:
                # concatenate umis
                umiseq = ''.join(umis)
            elif len(umis) == 1:
                umiseq = umis[0]
            # add umi to read name and write to outputfile
            writers[i].write(readnames[i] + ":" + umiseq + " " + namerests[i] + "\n")
            # remove umi and spacer from read seq. write remaining of read to outputfile
            if i == 0:
                # determine index for reads k <- 0 for r1, -1 for r2 or r3
                k = 0
            elif i > 0:
                k = -1
            
            # compute read length
            ReadLength[i].add(len(reads[k][1][UmiPositions[i] + UmiLength[i] + SpacerLength[i]:]))
            
            # write new fastqs
            writers[i].write(reads[k][1][UmiPositions[i] + UmiLength[i] + SpacerLength[i]:])
            writers[i].write(reads[k][2])
            writers[i].write(reads[k][3][UmiPositions[i] + UmiLength[i] + SpacerLength[i]:])
        
    # close all open files
    for i in writers:
        i.close()
    for i in fastqs:
        i.close()
    if spacer == True:
        if len(discarded) != 0:
            for i in discarded:
                i.close()
    
    print("Complete. Output written to {0}".format(outdir))
    
    # record read indo as json
    D = {'Total': Total, 'Correct': Correct, 'Incorrect': Incorrect}
    if len(ReadLength[0]) != 0:
        lr1 = list(ReadLength[0])[0]
        D['length_read1'] = lr1   
    if len(ReadLength[1]) != 0:
        lr2 = list(ReadLength[1])[0]
        D['length_read2'] = lr2
    return D, UmiSequences
