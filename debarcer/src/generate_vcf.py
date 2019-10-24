import sys
#import pysam
import configparser
import argparse
import operator
import functools
import csv
#from src.utilities import get_consecutive_items
import itertools

#from src.generate_consensus import ConsDataRow
#from src.generate_consensus import get_ref_seq
#from src.umi_error_correct import UMIGroup
#from src.generate_consensus import get_consensus_seq, get_uncollapsed_seq
#from src.handle_args import handle_arg
#from src.generate_consensus import temp_umi_table
#from src.generate_consensus import generate_uncollapsed
#from src.generate_consensus import generate_consensus
#from src.generate_consensus import memoize
#from src.generate_merge import check_consfile




#def write_rows(cons_data, f_size, contig, region_start, region_end, output_path, config):
#
#    with open("{}/{}:{}-{}.{}.vcf".format(output_path, contig, region_start, region_end, f_size), "w") as writer:
#
#        writer.write("##fileformat=VCFv4.2\n")
#        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
#        writer.write("##source=Debarcer2\n")
#        writer.write("##f_size={}\n".format(f_size))
#        
#        ## INFO/FILTER/FORMAT metadata
#        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
#        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
#        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
#        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
#        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
#        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
#        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
#        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")
#  
#        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
#        
#        for row in cons_data[f_size]:
#            writer.write(row)


#def vcf_output(cons_file, f_size, contig, region_start, region_end, output_path, config):
#    """Writes a .vcf consensus file."""
#        
#    with open("{}/{}:{}-{}.fsize{}.vcf".format(output_path, contig, region_start, region_end, f_size), "w") as writer:
#        
#        writer.write("##fileformat=VCFv4.2\n")
#        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
#        writer.write("##source=Debarcer2\n")
#        writer.write("##f_size={}\n".format(f_size))
#        
#        ## INFO/FILTER/FORMAT metadata
#        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
#        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
#        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
#        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
#        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
#        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
#        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
#        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")
#        
#        ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
#        all_threshold = float(config['REPORT']['percent_allele_threshold']) if config else 2.0
#        
#        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
#
#        for base_pos in range(region_start, region_end):
#            if base_pos in cons_data[f_size]:
#                
#                row = cons_data[f_size][base_pos]
#                ref = row.get_ref_info()
#                cons = row.get_cons_info()
#                stats = row.get_stats()
#
#                if stats['ref_freq'] <= ref_threshold:
#                        
#                    alleles = row.get_alleles(all_threshold)
#                    ref_bases = set([allele[0] for allele in alleles])
#                    ref_allele = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
#                    depths = row.impute_allele_depths()
#                    ref_depth = depths[ref_allele] if ref_allele in depths else 0
#                    alt_freqs = row.impute_allele_freqs(all_threshold)
#
#                    info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
#                        stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
#                    fmt_string = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency
#
#                    for ref_base in ref_bases:
#                        snips = []
#                        for allele in alleles:
#                            if allele[0] == ref_base:
#                                snips.append(allele)
#
#                        alt_string = ','.join( [allele[1] for allele in snips] )
#                        depth_string = ','.join( [str(depths[allele]) for allele in snips] )
#                        freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in snips] )
#                        smp_string = "{}:{}:{}".format(ref_depth, depth_string, freq_string)
#                        filt = "PASS" if any( [depths[alt] > 10 for alt in snips] ) else "a10"
#
#                        writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
#                            contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))
    


#def create_vcf_output(cons_file, f_size, cons_data, ref_seq, contig, region_start, region_end, output_path, config):
#    #base_positions, ref_bases = {}, {}
#    ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
#    all_threshold = float(config['REPORT']['percent_allele_threshold']) if config else 2.0
#
#	
#    #Read cons file
#    with open(cons_file, "r") as reader:
#        print("Reading cons file")
#
#        #f_size = 0
#        i = region_start
#        for line in reader:
#            ref_freq = line.split('\t')[13]
#
#            if i <= region_end:
#
#                #contig[i] = line.split('\t')[1] #### CHANGE: Added the following lines
#                base_positions[i] = line.split('\t')[2]
#                ref_bases[i] = line.split('\t')[3]
#                #rawdps[i] = line.split('\t')[10]
#                #consdps[i] = line.split('\t')[11]
#                #f_size[i] = line.split('\t')[12] ##### CHANGE: from [11] to [12]
# 
#                i=+1
#
#
#
#    with open("{}/{}:{}-{}.{}.vcf".format(output_path, contig, region_start, region_end, f_size), "w") as writer:
#        
#        writer.write("##fileformat=VCFv4.2\n")
#        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
#        writer.write("##source=Debarcer2\n")
#        writer.write("##f_size={}\n".format(f_size))
#        
#        ## INFO/FILTER/FORMAT metadata
#        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
#        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
#        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
#        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
#        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
#        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
#        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
#        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")
#        
#        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
#
#        for base_pos in range(region_start, region_end):
#            #if any( [base_pos in cons_data[f_size] for f_size in cons_data]):
#                #for f_size in cons_data:
#            if base_pos in cons_data[f_size]:
#
#                row = cons_data[f_size][base_pos]
#                ref = row.get_ref_info()
#                cons = row.get_cons_info()
#                stats = row.get_stats()
#
#
#                if stats['ref_freq'] <= ref_threshold:
#
#                    alleles = row.get_alleles(all_threshold)
#                    depths = row.impute_allele_depths()
#                    ref_bases = set([allele[0] for allele in alleles])
#                    ref_allele = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
#                    depths = row.impute_allele_depths()
#                    ref_depth = depths[ref_allele] if ref_allele in depths else 0
#                    alt_freqs = row.impute_allele_freqs(all_threshold)
#
#
#                    info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
#                        stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
#                    fmt_string = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency
#
#                    for ref_base in ref_bases:
#                        snips = []
#                        for allele in alleles:
#                            if allele[0] == ref_base:
#                                snips.append(allele)
#
#
#                        alt_string = ','.join( [allele[1] for allele in snips] )
#                        depth_string = ','.join( [str(depths[allele]) for allele in snips] )
#                        freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in snips] )
#                        smp_string = "{}:{}:{}".format(ref_depth, depth_string, freq_string)
#                        filt = "PASS" if any( [depths[alt] > 10 for alt in snips] ) else "a10"
#
#                        writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
#                            contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))



        







#def generate_vcf_output(cons_file, f_sizes, contig, region_start, region_end, output_path, config):
#    """(Main) generates VCF output file(s)."""
#
#    ## Get reference sequence for the region 
#    ref_seq = get_ref_seq(contig, region_start, region_end, config)
#
#    ## Parse table to extract VCF events
#    cons_data = parse_raw_table(cons_file, f_sizes)
#
#    ## Generate vcf files
#    for f_size in f_sizes:
#        if f_size in cons_data:
#            write_rows(cons_data, f_size, contig, region_start, region_end, output_path, config)




def check_consfile(cons_file):
    f = open(cons_file, "r")
    #reader = csv.reader(f, delimiter='\t')
    line = f.readline()
    if 'INTVL' in line:
        return True
    else:
        return False




def write_vcf(ID, config, contigs, f_sizes, base_positions, region_start, region_end, ref_base, alt_string, filt, info, fmt_string, smp_string, output_path, cons_is_merged):

    ref_threshold = 100.0
    #fam_idx = 0

    if not cons_is_merged:
        for i in range(len(contigs)):
            for size in f_sizes:
                if contigs[i][size] != None:
                    contig = contigs[i][size]
                    break

    #Create vcf file for a given base position for specified min fam size
    for f_size in f_sizes:

        if cons_is_merged:
            writer = open("{}/RUNID_{}.{}.vcf".format(output_path, ID, str(f_size)), "w")
        else:
            writer = open("{}/{}:{}-{}.{}.vcf".format(output_path, contig, region_start, region_end, str(f_size)), "w")

        writer.write("##fileformat=VCFv4.2\n")
        #writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
        
        writer.write("##reference={}\n".format('reference'))
                     
                     
        writer.write("##source=Debarcer2\n")
        writer.write("##f_size={}\n".format(f_size))

        ## INFO/FILTER/FORMAT metadata
        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth (ref allele depth, alt allele depth(s))\">\n")
        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Allele Frequency (ref allele frequency, alt allele freq(s))\">\n")
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for index in range(len(contigs)):
            if contigs[index][f_size] != None:
                writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                       contigs[index][f_size], base_positions[index][f_size], ".", ref_base[index][f_size], alt_string[index][f_size], "0", filt[index][f_size], info[index][f_size], fmt_string, smp_string[index][f_size]))
        #fam_idx += 1










######################################################
    
    
    
    
    
    
    
def vcf_output(cons_file, f_size, contig, region_start, region_end, output_path, config):
    """Writes a .vcf consensus file."""
        
    with open("{}/{}:{}-{}.fsize{}.vcf".format(output_path, contig, region_start, region_end, f_size), "w") as writer:
        
        writer.write("##fileformat=VCFv4.2\n")
        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
        writer.write("##source=Debarcer2\n")
        writer.write("##f_size={}\n".format(f_size))
        
        ## INFO/FILTER/FORMAT metadata
        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")
        
        ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
        all_threshold = float(config['REPORT']['percent_allele_threshold']) if config else 2.0
        
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for base_pos in range(region_start, region_end):
            if base_pos in cons_data[f_size]:
                
                row = cons_data[f_size][base_pos]
                ref = row.get_ref_info()
                cons = row.get_cons_info()
                stats = row.get_stats()

                if stats['ref_freq'] <= ref_threshold:
                        
                    alleles = row.get_alleles(all_threshold)
                    ref_bases = set([allele[0] for allele in alleles])
                    ref_allele = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
                    depths = row.impute_allele_depths()
                    ref_depth = depths[ref_allele] if ref_allele in depths else 0
                    alt_freqs = row.impute_allele_freqs(all_threshold)

                    info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
                        stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
                    fmt_string = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency

                    for ref_base in ref_bases:
                        snips = []
                        for allele in alleles:
                            if allele[0] == ref_base:
                                snips.append(allele)

                        alt_string = ','.join( [allele[1] for allele in snips] )
                        depth_string = ','.join( [str(depths[allele]) for allele in snips] )
                        freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in snips] )
                        smp_string = "{}:{}:{}".format(ref_depth, depth_string, freq_string)
                        filt = "PASS" if any( [depths[alt] > 10 for alt in snips] ) else "a10"

                        writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))
    
    
    
    
def parse_raw_table(cons_file, f_sizes):
    """Parses a .cons file generated by generate_consensus into VCF entries."""
    
    rows = {f_size: [] for f_size in f_sizes}

    with open(cons_file, "r") as reader:

        f_size = 0

        for line in reader:

            if line.startswith('#'):
                f_size = line.split('\t')[11]

            elif f_size in f_sizes:
                rows[f_size].append(line)

    return rows
    
    
    
def generate_vcf_output(cons_file, f_sizes, contig, region_start, region_end, output_path, config):
    """(Main) generates VCF output file(s)."""

    ## Get reference sequence for the region 
    ref_seq = get_ref_seq(contig, region_start, region_end, config)

    ## Parse table to extract VCF events
    cons_data = parse_raw_table(cons_file, f_sizes)

    ## Generate vcf files
    for f_size in f_sizes:
        if f_size in cons_data:
            write_rows(cons_data, f_size, contig, region_start, region_end, output_path, config)    
    
    
def parse_raw_table(cons_file, f_sizes):
    """Parses a .cons file generated by generate_consensus into VCF entries."""
    
    rows = {f_size: [] for f_size in f_sizes}

    with open(cons_file, "r") as reader:

        f_size = 0

        for line in reader:

            if line.startswith('#'):
                f_size = line.split('\t')[11]

            elif f_size in f_sizes:
                rows[f_size].append(line)

    return rows    
    


def get_consecutive_items(L):
    '''
    (list) -> generator
    
    :param L: List
    
    Return a generator with 1st and last item of consecutive items in L
    '''
    
    # remove duplicate and sort L
    L = sorted(set(L))
    for key, group in itertools.groupby(enumerate(L), lambda t: t[1] - t[0]):
        group = list(group)
        # return a generator with 1st and last item of consecutive items in L
        yield group[0][1], group[-1][1]    



def GetConsecutivePositions(consfile):
    '''
    (str) -> dict
    
    :param consfile: Path to the consensus file (merged or not)
    
    Returns a dictionary with list of consecutive positions (tuples with 1st and
    last positions included) per contig
    '''
    
    # create dict with all positions per contig
    positions = {}
    infile = open(consfile)
    header = infile.readline().rstrip().split('\t')
    # collect all positions per chromo
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            chromo, pos = line[0], int(line[1])
            if chromo not in positions:
                positions[chromo] = []
            positions[chromo].append(pos)
    infile.close()    
        
    # create a dict {contig: [(start1, end1), (start2, end2)]}
    conspos = {}
    for chromo in positions:
        conspos[chromo] = list(get_consecutive_items(positions[chromo]))
        
    return conspos

   
    
def GetConsData(consfile):
    '''
    (str) -> dict
    
    :param consfile: Path to the consensus file (merged or not)
    
    Returns a dictionary with consensus file info organized by chromo, region,
    and family size for each position
    '''

    # create a dict with list of consecutive positions per contig 
    # necessary for merged cons file
    conspos = GetConsecutivePositions(consfile) 

    # create a dict with consensus info per contig, region and umi family size
    # {contig :{(start, end): {fam: {pos: info}}}}
    data = {}

    infile = open(consfile)
    header = infile.readline().rstrip().split('\t')
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            contig = line[0]
            if contig not in data:
                data[contig] = {}
            # get the position
            pos = int(line[1])
            region = ''
            # find the region including pos
            for i in conspos[contig]:
                if pos in range(i[0], i[1] + 1):
                    region = i
                    break
            assert region != ''
            if region not in data[contig]:
                data[contig][region] = {}
            # get umi fam size              
            famsize = line[12]
            if famsize not in data[contig][region]:
                data[contig][region][famsize] = {}
            # collect info     
            assert pos not in data[contig][region][famsize]
            data[contig][region][famsize][pos] = line
    
    infile.close()
    return data
        
    
    
    
#############################################################    


# copy paste from generate_consensus    
    
#    row = cons_data[f_size][base_pos]
#                        ref = row.get_ref_info()
#                        cons = row.get_cons_info()
#                        stats = row.get_stats()
#
#                        if stats['ref_freq'] <= ref_threshold:
#
#                            alleles = row.get_alleles(all_threshold)
#                            ref_bases = set([allele[0] for allele in alleles])
#                            ref_allele = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
#                            depths = row.impute_allele_depths()
#                            ref_depth = depths[ref_allele] if ref_allele in depths else 0
#                            alt_freqs = row.impute_allele_freqs(all_threshold)
#
#                            info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
#                                stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
#                            fmt_string = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency
#
#                            for ref_base in ref_bases:
#
#                                #Handle error where ref_base is assigned to multiple ref_alleles
#                                if len(ref_base) > 1:
#                                    ref_base = ref_base[0]
#
#                                if (base_pos < 170837514) and (base_pos > 170837510):
#                                    print("--BASE_POS: "+str(base_pos)+" ref_base: "+ref_base+"--")
#
#                                snips = []
#                                for allele in alleles:
#                                    if allele[0] == ref_base:
#                                        snips.append(allele)
#
#                                alt_string = ','.join( [allele[1] for allele in snips] )
#                                depth_string = ','.join( [str(depths[allele]) for allele in snips] )
#                                freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in snips] )
#                                smp_string = "{}:{}:{}".format(ref_depth, depth_string, freq_string)
#                                filt = "PASS" if any( [depths[alt] > 10 for alt in snips] ) else "a10"
#
#                                #writer.write("# {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))
#
#

