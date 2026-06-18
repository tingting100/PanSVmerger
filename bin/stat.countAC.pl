#!/usr/bin/perl
use strict;
use warnings;

# USE: perl count_alleles.pl input.vcf.gz > allele_counts.txt
# OUT: CHROM  POS  REF  ALT_COUNT  TOTAL_ALLELES

my $vcf = $ARGV[0] or die "Usage: perl $0 <input.vcf|input.vcf.gz>\n";

my $fh;
if ($vcf =~ /\.gz$/) {
    open($fh, "-|", "gzip -dc $vcf") or die "Cannot open $vcf: $!";
} else {
    open($fh, "<", $vcf);
}

while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^#/;
    my @cols = split /\t/, $line;
    next if($cols[0] =~/chrM/ || $cols[0] =~/chrY/);
    my ($chrom, $pos, $ref, $alt_str, $info) = ($cols[0], $cols[1], $cols[3], $cols[4], $cols[7]);
    
    if ($info =~ /AC=([^;]+)/) {
        my $ac_str = $1;
        my @ac_vals = split /,/, $ac_str;
        my $alt_count = scalar @ac_vals;          
        my $total_alleles = $alt_count + 1;      
        
        print join("\t", $chrom, $pos, $ref, $alt_count, $total_alleles), "\n";
    } else {
        warn "No AC field at $chrom:$pos\n";
    }
}
close $fh;
