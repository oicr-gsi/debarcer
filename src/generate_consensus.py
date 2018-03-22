
import sys
import pysam
import configparser
import argparse
import operator
import functools
from src.get_consensus_seq import get_consensus_seq, get_uncollapsed_seq
from src.handle_args import handle_arg

class ConsDataRow:
    """Holds consensus data for one base position."""
    ref_info  = {}
    cons_info = {}
    stats     = {}
    
    def __init__(self, ref_info, cons_info, stats):
        self.ref_info  = ref_info
        self.cons_info = cons_info
        self.stats     = stats

    def impute_allele_depths(self):
        """Returns allele depths (dict of str tuples -> ints)."""
        return self.cons_info
    
    @functools.lru_cache(maxsize=2, typed=False)
    def impute_allele_freqs(self, threshold):
        """Returns allele frequencies (dict of str tuples -> floats)."""
        
        is_ref  = lambda allele: allele[0] is allele[1]
        
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


def generate_consensus(families, f_size, ref_seq, contig, region_start, region_end, bam_file, config):
    """Generates consensus data for the given family size and region."""

    ## Keys: each base position in the region
    ## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
    consensus_seq = get_consensus_seq(families, ref_seq, contig, region_start, region_end, bam_file, config)

    percent_threshold = float(config['SETTINGS']['percent_consensus_threshold']) if config else 70.0
    count_threshold   = int(config['SETTINGS']['count_consensus_threshold']) if config else 1

    cons_data = {}

    for base_pos in range(region_start, region_end):

        ref_base = ref_seq[base_pos-region_start]

        if base_pos in consensus_seq:

            consensuses = {}
            raw_depth   = 0
            min_fam     = max([sum(consensus_seq[base_pos][fam].values()) 
                            for fam in consensus_seq[base_pos]]) 
            
            for family in consensus_seq[base_pos]:
                        
                cons_allele  = max(consensus_seq[base_pos][family].items(), key = operator.itemgetter(1))[0]
                cons_denom   = sum(consensus_seq[base_pos][family].values())
                cons_percent = (consensus_seq[base_pos][family][cons_allele]/cons_denom) * 100
                
                raw_depth += cons_denom
    
                if cons_percent >= percent_threshold and consensus_seq[base_pos][family][cons_allele] >= count_threshold:
                    
                    if cons_allele in consensuses:
                        consensuses[cons_allele] += 1

                    else:
                        consensuses[cons_allele] = 1
                        
                    if sum(consensus_seq[base_pos][family].values()) < min_fam:
                        min_fam = sum(consensus_seq[base_pos][family].values())
            
            cons_depth = len(consensus_seq[base_pos])
            mean_fam   = sum( [sum(consensus_seq[base_pos][fam].values()) 
                            for fam in consensus_seq[base_pos]] ) / len(consensus_seq[base_pos])
            ref_freq   = (consensuses[(ref_base, ref_base)] / cons_depth) * 100 if (ref_base, ref_base) in consensuses else 0
            
            ref_info  = {"contig": contig, "base_pos": base_pos, "ref_base": ref_base}
            cons_info = consensuses
            stats     = {"rawdp": raw_depth, "consdp": cons_depth, "min_fam": min_fam, "mean_fam": mean_fam, "ref_freq": ref_freq}
                    
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
                
            depth    = sum(uncollapsed_seq[base_pos].values())
            ref_freq = (uncollapsed_seq[base_pos][ref_base] / depth) * 100 if ref_base in uncollapsed_seq[base_pos] else 0
            get_cons = lambda base: uncollapsed_seq[base_pos][base] if base in uncollapsed_seq[base_pos] else 0    
            
            ref_info  = {"contig": contig, "base_pos": base_pos, "ref_base": ref_base}
            cons_info = uncollapsed_seq[base_pos]
            stats     = {"rawdp": depth, "consdp": None, "min_fam": 0, "mean_fam": 0, "ref_freq": ref_freq}
            
            row = ConsDataRow(ref_info, cons_info, stats)
            cons_data[base_pos] = row
    
    return cons_data
    

def tabular_output(cons_data, contig, region_start, region_end, output_path, config):
    """Writes a .txt consensus file with tabular columns."""
    
    with open("{}/{}:{}-{}.cons".format(output_path, contig, region_start, region_end), "w") as writer:
        
        for base_pos in range(region_start, region_end):
            for f_size in cons_data:                
                if base_pos in cons_data[f_size]:
                    
                    ref   = cons_data[f_size][base_pos].get_ref_info()
                    cons  = cons_data[f_size][base_pos].get_cons_info()
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
                    
                    writer.write("{}\t{}\t{}\t".format(
                        ref['contig'], ref['base_pos'], ref['ref_base']))
                    writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                        counts['A'], counts['C'], counts['G'], counts['T'], counts['I'], counts['D'], counts['N']))
                    writer.write("{}\t{}\t{}\t{}\n".format(
                        stats['rawdp'], stats['consdp'], stats['min_fam'], stats['ref_freq']))
                    

def vcf_output(cons_data, f_size, ref_seq, contig, region_start, region_end, output_path, config):
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
        writer.write("##FILTER=<ID=a10,Number=0,Type=Flag,Description=\"Alt allele depth below 10\">\n")
        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")

        writer.write("#CHROM\tPOS\t\tID\tREF\tALT\t\tQUAL\tFILTER\tINFO\t\t\tFORMAT\t\tSAMPLE\n") ## TODO multiple samples?
        
        ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
        all_threshold = float(config['REPORT']['percent_allele_threshold']) if config else 2.0
        
        for base_pos in range(region_start, region_end):
                if base_pos in cons_data[f_size]:
                    
                    row   = cons_data[f_size][base_pos]
                    ref   = row.get_ref_info()
                    cons  = row.get_cons_info()
                    stats = row.get_stats()
                    
                    if stats['ref_freq'] <= ref_threshold:
                        
                        alleles     = row.get_alleles(all_threshold)
                        ref_allele  = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
                        ref_string  = ','.join( [allele[0] for allele in alleles] )
                        alt_string  = ','.join( [allele[1] for allele in alleles] )
                        depths      = row.impute_allele_depths()
                        ref_depth   = depths[ref_allele] if ref_allele in depths else 0
                        alt_depths  = ','.join( [str(depths[allele]) for allele in alleles] )
                        alt_freqs   = row.impute_allele_freqs(all_threshold)
                        freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in alt_freqs] )
                    
                        filt = "PASS" if any( [depths[alt] > 10 for alt in alleles] ) else "a10"
                        info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
                            stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
                        fmt  = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency
                        smp  = "{}:{}:{}".format(ref_depth, alt_depths, freq_string)
                        
                        writer.write("{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\n".format(
                            contig, base_pos, None, ref_string, alt_string, None, filt, info, fmt, smp))
        

def generate_consensus_output(contig, region_start, region_end, bam_file, tally_file, output_path, config):
    """(Main) generates tabular and VCF consensus output files."""

    ## Get reference sequence
    with pysam.FastaFile(config['PATHS']['reference_file']) as reader:
        ref_seq = reader.fetch(contig, region_start, region_end).upper()
        
    ## Lists of UMI+Pos pairs with count >= f_size
    families = {}
    f_sizes  = [int(n) for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else [1, 2, 5, 10]

    with open(tally_file, "r") as reader:
        lines = reader.readlines()
    
    for line in lines[1:]:
        umi, pos, count, *rest = line.split("\t")
    
        for f_size in f_sizes:
        
            if(int(count) >= f_size):
            
                if f_size not in families:
                    families[f_size] = {}
                
                families[f_size][umi + pos] = 0

    ## Get consensus data for each f_size + uncollapsed data
    cons_data    = {}
    cons_data[0] = generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, config)

    for f_size in families:
        cons_data[f_size] = generate_consensus(
            families[f_size], f_size, ref_seq, contig, region_start, region_end, bam_file, config)
    
    ## Output
    tabular_output(cons_data, contig, region_start, region_end, output_path, config)
    
    for f_size in families:
        vcf_output(cons_data, f_size, ref_seq, contig, region_start, region_end, output_path, config)


if __name__=="__main__":
    
    ## Argument + config parsing and error handling
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam_file',    help='Path to your BAM file.')
    parser.add_argument('-r', '--region',      help='Region to analyze (string of the form chrX:posA-posB).')
    parser.add_argument('-o', '--output_path', help='Path to write output files to.')
    parser.add_argument('-c', '--config',      help='Path to your config file.')
    parser.add_argument('-t', '--tally',       help='Path to your tally (output of UMI_count.py).')

    args = parser.parse_args()
    config_file = args.config

    if config_file:
        config = configparser.ConfigParser()
        config.read(config_file)
    else:
        config = None
    
    region = args.region
    if any(x not in region for x in ["chr", ":", "-"]):
        raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
        sys.exit(1)

    contig       = region.split(":")[0]
    region_start = int(region.split(":")[1].split("-")[0])
    region_end   = int(region.split(":")[1].split("-")[1])

    bam_file    = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None, 
                    'No BAM file provided in args or config.')
    output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
                    'No output path provided in args or config.')
    tally_file  = handle_arg(args.tally, output_path + '/' + region + '.tally' if config else None, 
                    'No tally file provided.')

    ## Output
    generate_consensus_output(contig, region_start, region_end, bam_file, tally_file, output_path, config)
    
    
        
    
        
    
            
            
            
    
    