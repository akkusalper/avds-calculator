#!/usr/bin/env python3
"""
Create synthetic test data for AVDS testing
Generates a small BAM file and VCF file for testing purposes
"""

import pysam
import os

def create_test_bam(output_path, reference_name="chr22", reference_length=1000000):
    """
    Create a simple test BAM file with synthetic reads.

    Creates reads covering a few variant positions with controlled quality metrics.
    """
    # Define variant positions and their properties
    variants = [
        # (pos, ref, alt, vaf, depth, strand_balance, quality)
        (100000, 'A', 'G', 0.50, 100, 0.50, 40),  # Good heterozygous variant
        (200000, 'C', 'T', 0.30, 80, 0.45, 35),   # Lower VAF variant
        (300000, 'G', 'A', 0.50, 120, 0.10, 38),  # Strand bias artifact
        (400000, 'T', 'C', 0.50, 50, 0.50, 25),   # Lower quality
        (500000, 'A', 'T', 0.05, 200, 0.50, 40),  # Very low VAF
    ]

    # Create header
    header = {
        'HD': {'VN': '1.0', 'SO': 'coordinate'},
        'SQ': [{'LN': reference_length, 'SN': reference_name}],
        'RG': [{'ID': 'test', 'SM': 'test_sample', 'PL': 'ILLUMINA'}]
    }

    # Open output BAM file
    with pysam.AlignmentFile(output_path, "wb", header=header) as outf:
        read_id = 0

        for pos, ref, alt, vaf, depth, strand_balance, base_qual in variants:
            # Calculate number of reads for each allele
            n_alt = int(depth * vaf)
            n_ref = depth - n_alt

            # Calculate strand distribution
            n_forward = int(depth * strand_balance)
            n_reverse = depth - n_forward

            # Distribute alt reads between strands
            n_alt_forward = int(n_alt * strand_balance)
            n_alt_reverse = n_alt - n_alt_forward

            # Create reference reads
            for i in range(n_ref):
                read_id += 1
                is_forward = i < (n_ref * strand_balance)
                read = create_read(
                    read_id, reference_name, pos, ref,
                    is_forward, base_qual, is_alt=False, ref_base=ref
                )
                outf.write(read)

            # Create alternative reads
            for i in range(n_alt):
                read_id += 1
                is_forward = i < n_alt_forward
                read = create_read(
                    read_id, reference_name, pos, alt,
                    is_forward, base_qual, is_alt=True, ref_base=ref
                )
                outf.write(read)

    # Create index
    pysam.index(output_path)
    print(f"Created BAM file: {output_path}")
    print(f"Created BAM index: {output_path}.bai")


def create_read(read_id, chrom, var_pos, base, is_forward, base_qual, is_alt=True, ref_base='A'):
    """Create a synthetic aligned read"""
    read_length = 150

    # Position read so variant is at position 74 (0-indexed)
    # var_pos is 1-based, so var_pos-1 is 0-based
    # We want: read_start + 74 = var_pos - 1
    # Therefore: read_start = var_pos - 1 - 74 = var_pos - 75
    read_start = var_pos - 1 - 74
    read_end = read_start + read_length

    # Create sequence with variant at position 74 (0-indexed)
    # This ensures reads are properly matched to ref/alt
    sequence = ref_base * 74 + base + ref_base * 75

    # Create quality scores
    qualities = [base_qual] * read_length

    # Create read
    read = pysam.AlignedSegment()
    read.query_name = f"read_{read_id}"
    read.query_sequence = sequence
    read.flag = 0 if is_forward else 16  # 16 = reverse strand
    read.reference_id = 0  # First (and only) chromosome
    read.reference_start = read_start
    read.mapping_quality = 60
    read.cigar = [(0, read_length)]  # 150M (match)
    read.query_qualities = qualities
    read.set_tag('AS', read_length)  # Alignment score
    read.set_tag('NM', 0 if is_alt else 0)  # Edit distance
    read.set_tag('RG', 'test')

    return read


def create_test_vcf(output_path, reference_name="chr22"):
    """Create a simple test VCF file with the same variants as the BAM"""
    variants = [
        (100000, 'A', 'G'),
        (200000, 'C', 'T'),
        (300000, 'G', 'A'),
        (400000, 'T', 'C'),
        (500000, 'A', 'T'),
    ]

    # Create header
    header = pysam.VariantHeader()
    header.add_line(f'##contig=<ID={reference_name},length=1000000>')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_sample('test_sample')

    # Write VCF
    with pysam.VariantFile(output_path, 'w', header=header) as vcf:
        for pos, ref, alt in variants:
            record = vcf.new_record()
            record.chrom = reference_name
            record.pos = pos
            record.id = '.'
            record.ref = ref
            record.alts = (alt,)
            record.qual = None
            record.filter.add('PASS')
            record.samples['test_sample']['GT'] = (0, 1)  # Heterozygous
            vcf.write(record)

    print(f"Created VCF file: {output_path}")


def main():
    """Main function to create test data"""
    # Output directory
    output_dir = os.path.dirname(os.path.abspath(__file__))

    # Create test files
    bam_path = os.path.join(output_dir, "test.bam")
    vcf_path = os.path.join(output_dir, "test.vcf.gz")

    print("Creating synthetic test data...")
    print(f"Output directory: {output_dir}")
    print()

    # Create BAM
    create_test_bam(bam_path)
    print()

    # Create VCF
    create_test_vcf(vcf_path)

    # Compress and index VCF
    pysam.tabix_index(vcf_path, preset='vcf', force=True)
    print(f"Created VCF index: {vcf_path}.tbi")
    print()

    print("=" * 80)
    print("Test data created successfully!")
    print("=" * 80)
    print(f"BAM file: {bam_path}")
    print(f"VCF file: {vcf_path}")
    print()
    print("You can now test AVDS with:")
    print(f"python3 ../avds_pipeline.py -v {vcf_path} -b {bam_path} -o avds_results.tsv")
    print("=" * 80)


if __name__ == '__main__':
    main()
