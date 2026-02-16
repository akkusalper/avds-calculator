#!/usr/bin/env python3
"""Debug test data creation"""

# Simulate first variant
pos, ref, alt, vaf, depth, strand_balance, base_qual = (100000, 'A', 'G', 0.50, 100, 0.50, 40)

# Calculate number of reads for each allele
n_alt = int(depth * vaf)
n_ref = depth - n_alt

# Calculate strand distribution
n_forward = int(depth * strand_balance)
n_reverse = depth - n_forward

# Distribute alt reads between strands
n_alt_forward = int(n_alt * strand_balance)
n_alt_reverse = n_alt - n_alt_forward

print(f"Variant: {ref}>{alt} at position {pos}")
print(f"VAF: {vaf}, Depth: {depth}")
print(f"\nCalculated:")
print(f"  n_alt: {n_alt}")
print(f"  n_ref: {n_ref}")
print(f"  n_forward (total): {n_forward}")
print(f"  n_reverse (total): {n_reverse}")
print(f"  n_alt_forward: {n_alt_forward}")
print(f"  n_alt_reverse: {n_alt_reverse}")

print(f"\nExpected alt reads to create: {n_alt}")
print(f"Expected ref reads to create: {n_ref}")
