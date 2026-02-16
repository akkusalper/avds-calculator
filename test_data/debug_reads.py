#!/usr/bin/env python3
import pysam

bam = pysam.AlignmentFile('test.bam', 'rb')

print('Checking variant at chr22:100000 (A>G)')
print('='*60)

n_total = 0
n_alt_g = 0
n_ref_a = 0
n_other = 0

for read in bam.fetch('chr22', 99999, 100001):
    n_total += 1
    try:
        # Get position in read
        pos_0based = 100000 - 1
        ref_pos = read.reference_start
        query_idx = 0

        query_pos = None
        for op, length in read.cigartuples:
            if op == 0:  # M
                if ref_pos <= pos_0based < ref_pos + length:
                    query_pos = query_idx + (pos_0based - ref_pos)
                    break
                ref_pos += length
                query_idx += length

        if query_pos is not None:
            base = read.query_sequence[query_pos]
            if base == 'G':
                n_alt_g += 1
            elif base == 'A':
                n_ref_a += 1
            else:
                n_other += 1

            if n_total <= 5:
                print(f"Read {n_total}: base={base}, forward={not read.is_reverse}, mapq={read.mapping_quality}")
    except Exception as e:
        print(f"Error processing read: {e}")

print(f'\nTotal reads: {n_total}')
print(f'Ref (A): {n_ref_a}')
print(f'Alt (G): {n_alt_g}')
print(f'Other: {n_other}')

bam.close()
