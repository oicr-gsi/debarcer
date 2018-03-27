
import pysam
import configparser


def add_base(mode, seq, pos, family, allele):
    """(Helper) Adds a base to a seq dictionary."""

    added = False

    while not added and mode == "uncollapsed":
        if pos in seq:
            if allele in seq[pos]:
                seq[pos][allele] += 1

            else:
                seq[pos][allele] = 1

            added = True

        else:
            seq[pos] = {}

    while not added and mode == "consensus":
        if pos in seq:
            if family in seq[pos]:
                if allele in seq[pos][family]:
                    seq[pos][family][allele] += 1

                else:
                    seq[pos][family][allele] = 1

                added = True

            else:
                seq[pos][family] = {}

        else:
            seq[pos] = {}


def get_consensus_seq(families, ref_seq, contig, region_start, region_end, bam_file, config):
    """
    Returns a nested dictionary representing counts of each base in each family at each base pos'n.
     - Keys: each base position in the region
     - Values: tables of A,C,G,T (etc) counts from each UMI+Pos family
    """
    consensus_seq = {}

    pos_threshold = int(config['SETTINGS']['umi_family_pos_threshold']) if config else 10

    # Need to reconsider this for pileup, significant performance improvements possible
    ## read_threshold = int(config['SETTINGS']['max_family_size']) if config else 50

    with pysam.AlignmentFile(bam_file, "rb") as reader:

        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=1000000):

            pos = pileupcolumn.reference_pos
            if pos >= region_start and pos < region_end:

                for read in pileupcolumn.pileups:

                    read_data = read.alignment
                    read_name = read_data.query_name
                    start_pos = read_data.reference_start
                    umi = read_name.split(":")[-1]

                    for family in families:

                        if umi in family[0] and abs(family[1] - start_pos) <= pos_threshold:

                            family_key = family[0][0] + str(family[1]) ## unique family identifier

                            ref_pos = pos - region_start

                            if not read.is_del and not read.indel:
                                ref_base = ref_seq[ref_pos]
                                alt_base = read_data.query_sequence[read.query_position]

                            # Next position is an insert (current base is ref)
                            elif read.indel > 0:
                                ref_base = ref_seq[ref_pos]
                                alt_base = read_data.query_sequence[
                                    read.query_position:read.query_position + abs(read.indel) + 1]

                            # Next position is a deletion (current base + next bases are ref)
                            elif read.indel < 0:
                                ref_base = ref_seq[ref_pos:ref_pos + abs(read.indel) + 1]
                                alt_base = read_data.query_sequence[read.query_position]

                            add_base(mode="consensus", seq=consensus_seq, pos=pos,
                                     family=family_key, allele=(ref_base, alt_base))

                            break

    return consensus_seq


def get_uncollapsed_seq(ref_seq, contig, region_start, region_end, bam_file, config):
    """
    Returns a nested dictionary representing counts of each base at each base pos'n.
     - Keys: each base position in the region
     - Values: tables of A,T,C,G (etc) counts
    """

    uncollapsed_seq = {}

    with pysam.AlignmentFile(bam_file, "rb") as reader:

        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=1000000):

            pos = pileupcolumn.reference_pos
            if pos >= region_start and pos < region_end:

                for read in pileupcolumn.pileups:

                    if not read.is_del and not read.indel:
                        ref_base = ref_seq[pos - region_start]
                        alt_base = read.alignment.query_sequence[read.query_position]

                    # Next position is an insert (current base is ref)
                    elif read.indel > 0:
                        ref_base = ref_seq[pos - region_start]
                        alt_base = read.alignment.query_sequence[
                            read.query_position:read.query_position + abs(read.indel) + 1]

                    # Next position is a deletion (current base + next bases are ref)
                    elif read.indel < 0:
                        ref_base = ref_seq[read.query_position:read.query_position + abs(read.indel) + 1]
                        alt_base = read.alignment.query_sequence[read.query_position]

                    add_base(mode="uncollapsed", seq=uncollapsed_seq, pos=pos, family=None, allele=(ref_base, alt_base))

    return uncollapsed_seq
