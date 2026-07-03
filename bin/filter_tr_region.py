#!/usr/bin/env python3
"""
filter_tr_region.py - Filter VCF variants that fall in BED regions

Usage:
    python filter_tr_region.py --vcf input.vcf --bed regions.bed --out filtered.vcf
    python filter_tr_region.py --vcf input.vcf.gz --bed regions.bed --out filtered.vcf.gz
"""

import argparse
import gzip
import sys

def open_file(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)

def read_bed_regions(bed_file):
    regions = []
    with open_file(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                regions.append((chrom, start, end))
    return regions

def is_in_bed(chrom, pos, regions)
    for r_chrom, r_start, r_end in regions:
        if chrom == r_chrom and r_start <= pos <= r_end:
            return True
    return False

def main():
    parser = argparse.ArgumentParser(description='Filter VCF variants in BED regions')
    parser.add_argument('--vcf', required=True, help='Input VCF file (supports .gz)')
    parser.add_argument('--bed', required=True, help='BED file with regions to filter out')
    parser.add_argument('--out', required=True, help='Output VCF file (supports .gz)')
    args = parser.parse_args()


    print(f"Reading BED regions from {args.bed}...", file=sys.stderr)
    regions = read_bed_regions(args.bed)
    print(f"Loaded {len(regions)} regions", file=sys.stderr)

   
    vcf_in = open_file(args.vcf, 'r')
    vcf_out = open_file(args.out, 'w')

   
    total = 0
    filtered = 0

    for line in vcf_in:
        # keep header
        if line.startswith('#'):
            vcf_out.write(line)
            continue
        
       
        parts = line.split('\t')
        if len(parts) < 2:
            vcf_out.write(line)
            continue
        
        chrom = parts[0]
        pos = int(parts[1])
        total += 1
        
        
        if is_in_bed(chrom, pos, regions):
            filtered += 1
        else:
            vcf_out.write(line)

    vcf_in.close()
    vcf_out.close()

    kept = total - filtered
    print(f"\nSummary:", file=sys.stderr)
    print(f"  Total variants: {total:,}", file=sys.stderr)
    print(f"  Filtered out:   {filtered:,} ({filtered/total*100:.1f}%)", file=sys.stderr)
    print(f"  Kept:           {kept:,} ({kept/total*100:.1f}%)", file=sys.stderr)
    print(f"\nOutput written to: {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
