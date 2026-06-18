#!/usr/bin/env python3
"""
Extract all AF values from VCF for density plotting.
Usage: python extract_af.py --input_vcf input.vcf.gz --output_af af_values.txt
"""

import gzip
import argparse

def extract_af(vcf_file, output_file):
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    
    af_list = []
    
    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            info = fields[7]
            for item in info.split(';'):
                if item.startswith('AF='):
                    af_str = item.split('=')[1]
                    for af in af_str.split(','):
                        try:
                            af_list.append(float(af))
                        except ValueError:
                            continue
                    break
    
    with open(output_file, 'w') as f:
        f.write("AF\n")
        for af in af_list:
            f.write(f"{af}\n")
    
    print(f"Extracted {len(af_list)} AF values to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf", required=True, help="Input VCF (.vcf or .vcf.gz)")
    parser.add_argument("--output_af", required=True, help="Output file with AF values")
    args = parser.parse_args()
    extract_af(args.input_vcf, args.output_af)
