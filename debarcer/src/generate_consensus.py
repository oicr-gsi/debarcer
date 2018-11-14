
import sys
import pysam
import configparser
import argparse
import operator
import functools
from src.umi_error_correct import UMIGroup
from src.get_consensus_seq import get_consensus_seq, get_uncollapsed_seq
from src.handle_args import handle_arg
from src.get_ref_seq import get_ref_seq

class ConsDataRow:
    """Holds consensus data for one base position."""
    ref_info  = {}
    cons_info = {}
    stats     = {}
    
    def __init__(self, ref_info, cons_info, stats):
        self.ref_info = ref_info
        self.cons_info = cons_info
        self.stats = stats

    def impute_allele_depths(self):
        """Returns allele depths (dict of str tuples -> ints)."""
        return self.cons_info
    
    @functools.lru_cache(maxsize=2, typed=False)
    def impute_allele_freqs(self, threshold):
        """Returns allele frequencies (dict of str tuples -> floats)."""
        
        is_ref = lambda allele: allele[0] is allele[1]
        
        freqs = {}
        for allele in self.cons_info:
            freq = (self.cons_info[allele] / sum(self.cons_info.values())) * 100
            
            if not is_ref(allele) and freq > threshold:
                freqs[allele] = freq
        
        return freqs
    
    def get_alleles(self, threshold):
        """Returns alt alleles with their associated refs (list of str tuples)."""
        
        freqs   = self.impute_allele_freqs(threshold)
        alleles = []
        
        for allele in self.cons_info:
            if allele in freqs:
                alleles.append(allele)
        
        return alleles

    def get_ref_info(self):
        return self.ref_info
    
    def get_cons_info(self):
        return self.cons_info

    def get_stats(self):
        return self.stats


def generate_consensus(umi_table, f_size, ref_seq, contig, region_start, region_end, bam_file, config):
    """Generates consensus data for the given family size and region."""

    ## Keys: each base position in the region
    ## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
    consensus_seq = get_consensus_seq(umi_table, f_size, ref_seq, contig, region_start, region_end, bam_file, config)

    percent_threshold = float(config['SETTINGS']['percent_consensus_threshold']) if config else 70.0
    count_threshold = int(config['SETTINGS']['count_consensus_threshold']) if config else 1

    cons_data = {}

    for base_pos in range(region_start, region_end):

        ref_base = ref_seq[base_pos-region_start]

        if base_pos in consensus_seq:

            consensuses = {}
            raw_depth = 0
            min_fam = min([sum(consensus_seq[base_pos][fam].values()) for fam in consensus_seq[base_pos]]) 
            
            for family in consensus_seq[base_pos]:
                        
                cons_allele = max(consensus_seq[base_pos][family].items(), key = operator.itemgetter(1))[0]
                cons_denom = sum(consensus_seq[base_pos][family].values())
                cons_percent = (consensus_seq[base_pos][family][cons_allele]/cons_denom) * 100
                
                raw_depth += cons_denom
    
                if cons_percent >= percent_threshold and consensus_seq[base_pos][family][cons_allele] >= count_threshold:
                    
                    if cons_allele in consensuses:
                        consensuses[cons_allele] += 1

                    else:
                        consensuses[cons_allele] = 1
            
            cons_depth = sum(consensuses.values())
            mean_fam = sum( [sum(consensus_seq[base_pos][fam].values()) 
                            for fam in consensus_seq[base_pos]] ) / len(consensus_seq[base_pos])
            ref_freq = (consensuses[(ref_base, ref_base)] / cons_depth) * 100 if (ref_base, ref_base) in consensuses else 0
            
            ref_info = {"contig": contig, "base_pos": base_pos + 1, "ref_base": ref_base}
            cons_info = consensuses
            stats = {"rawdp": raw_depth, "consdp": cons_depth, "min_fam": min_fam, "mean_fam": mean_fam, "ref_freq": ref_freq}
                    
            row = ConsDataRow(ref_info, cons_info, stats)
            cons_data[base_pos] = row
                    

    return cons_data


def generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, config):
    """Generates uncollapsed consensus data for the given family size and region."""
    
    ## Keys: each base position in the region
    ## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
    uncollapsed_seq = get_uncollapsed_seq(ref_seq, contig, region_start, region_end, bam_file, config)
    
    cons_data = {}
    
    for base_pos in range(region_start, region_end):
            
        ref_base = ref_seq[base_pos-region_start]
            
        if base_pos in uncollapsed_seq:
                
            depth = sum(uncollapsed_seq[base_pos].values())
            
            if (ref_base, ref_base) in uncollapsed_seq[base_pos]:
                ref_freq = (uncollapsed_seq[base_pos][(ref_base, ref_base)] / depth) * 100
            else:
                ref_freq = 0

            ref_info = {"contig": contig, "base_pos": base_pos, "ref_base": ref_base}
            cons_info = uncollapsed_seq[base_pos]
            stats = {"rawdp": depth, "consdp": depth, "min_fam": 0, "mean_fam": 0, "ref_freq": ref_freq}
            
            row = ConsDataRow(ref_info, cons_info, stats)
            cons_data[base_pos] = row
    
    return cons_data


def raw_table_output(cons_data, ref_seq, contig, region_start, region_end, output_path, config):
    """Writes a long-form consensus file for every event detected in the collapsed data."""

    ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
    all_threshold = float(config['REPORT']['percent_allele_threshold']) if config else 2.0

    with open("{}/{}:{}-{}.cons".format(output_path, contig, region_start, region_end), "w") as writer:

        writer.write("CHROM\tPOS\tREF\t'A's\t'C's\t'G's\t'T's\t'I's\t'D's\t'N's\tRAWDP\tCONSDP\tFAM\tREF_FREQ\tMEAN_FAM\n") ##Header
        
        for base_pos in range(region_start, region_end):

            ref_base = ref_seq[base_pos-region_start]
            if any( [base_pos in cons_data[f_size] for f_size in cons_data] ):

                for f_size in cons_data:
                    if base_pos in cons_data[f_size]:

                        ref = cons_data[f_size][base_pos].get_ref_info()
                        cons = cons_data[f_size][base_pos].get_cons_info()
                        stats = cons_data[f_size][base_pos].get_stats()

                        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'I': 0, 'D': 0, 'N': 0}
                        for allele in cons:
                            # ref > 1 => deletion
                            if len(allele[0]) > 1:
                                counts['D'] += cons[allele]
                                
                            # allele > 1 => insertion
                            elif len(allele[1]) > 1:
                                counts['I'] += cons[allele]
                                
                            else:
                                counts[allele[1]] += cons[allele]


                        writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            contig, base_pos, ref_base, counts['A'], counts['C'], counts['G'], counts['T'], 
                            counts['I'], counts['D'], counts['N'],
                            stats['rawdp'], stats['consdp'], f_size, stats['ref_freq'], stats['mean_fam']))

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

                                #writer.write("# {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))


def temp_umi_table(contig, region_start, region_end, bam_file, config):
    """Makes a UMI table with no error correction (only when one is not provided)."""

    umi_table = {}
    pos_threshold = int(config['SETTINGS']['umi_family_pos_threshold']) if config else 10

    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:

        for read in bam_reader.fetch(contig, region_start, region_end):

            umi = read.query_name.split(':')[-1]
            pos = read.reference_start

            if umi in umi_table:

                umi_group = umi_table[umi]
                families = umi_group.families

                if pos in families:
                    umi_group.add(pos)

                else:
                    closest = umi_table[umi].getClosest(pos, pos_threshold)
                    if closest:
                        umi_group.add(closest)
                    else:
                        umi_group.addNew(pos)

            else:
                umi_table[umi] = UMIGroup(key=umi)
                umi_table[umi].addNew(pos)

    return umi_table



def memoize(func):
    cache = dict()

    def memoized_func(*args):
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result

    return memoized_func



def generate_consensus_output(contig, region_start, region_end, bam_file, umi_table, output_path, config):
    """(Main) generates consensus output file."""

    ## Make a stand-in umi_table if one is not provided (no error correction)
    if not umi_table:
        print("Building temporary UMI table...")
        umi_table = temp_umi_table(contig, region_start, region_end, bam_file, config)
        
    ## Lists of umi families with count >= f_size
    f_sizes = [int(n) for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else [1, 2, 5]

    ## Get reference sequence for the region 
    print("Getting reference sequence...")
    ref_seq = get_ref_seq(contig, region_start, region_end, config)

    ## Get consensus data for each f_size + uncollapsed data
    print("Building consensus data...")
    cons_data = {}
    
    """
    memoized_uncollapsed = memoize(generate_uncollapsed)
    cons_data[0] = memoized_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, config)

    memoized_collapsed = memoize(generate_consensus)


    for f_size in f_sizes:
        cons_data[f_size] = memoized_collapsed(umi_table, f_size, ref_seq, contig, region_start, region_end, bam_file, config)
    """

    cons_data[0] = generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, config)

    for f_size in f_sizes:
        cons_data[f_size] = generate_consensus(
            umi_table, f_size, ref_seq, contig, region_start, region_end, bam_file, config)
  

    ## Output
    print("Writing output...")
    raw_table_output(cons_data, ref_seq, contig, region_start, region_end, output_path, config)
    

if __name__=="__main__":
    
    ## Argument + config parsing and error handling
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam_file', help='Path to your BAM file.')
    parser.add_argument('-r', '--region', help='Region to analyze (string of the form chrX:posA-posB).')
    parser.add_argument('-o', '--output_path', help='Path to write output files to.')
    parser.add_argument('-c', '--config', help='Path to your config file.')
    parser.add_argument('-t', '--tally', help='Path to your tally (output of UMI_count.py).')

    args = parser.parse_args()

    if args.config:
        config = configparser.ConfigParser()
        config.read(config_file)
    else:
        config = None
    
    region = args.region
    if any(x not in region for x in ["chr", ":", "-"]):
        raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
        sys.exit(1)

    contig = region.split(":")[0]
    region_start = int(region.split(":")[1].split("-")[0])
    region_end = int(region.split(":")[1].split("-")[1])

    bam_file = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None, 
                    'No BAM file provided in args or config.')
    output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
                    'No output path provided in args or config.')
    tally_file = handle_arg(args.tally, output_path + '/' + region + '.tally' if config else None, 
                    'No tally file provided.')

    ## Output
    generate_consensus_output(contig, region_start, region_end, bam_file, tally_file, output_path, config)
