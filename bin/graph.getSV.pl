#!/usr/bin/perl -w
use strict;
use warnings;

my $vcf_file = $ARGV[0];
my $out_file = $ARGV[1];

die "Usage: perl $0 <input.vcf> <output.vcf>\n" unless $vcf_file && $out_file;

my $in_fh;
if ($vcf_file =~ /\.gz$/) {
    open $in_fh, "gzip -cd $vcf_file |" or die "Cannot open $vcf_file: $!";
} else {
    open $in_fh, "<", $vcf_file or die "Cannot open $vcf_file: $!";
}

my $out_fh;
if ($out_file =~ /\.gz$/) {
    open $out_fh, "| gzip -c > $out_file" or die "Cannot write to $out_file: $!";
} else {
    open $out_fh, ">", $out_file or die "Cannot write to $out_file: $!";
}

my $total = 0;
my $sv_count = 0;

while (<$in_fh>) {
    chomp;
    
    if (/^#/) {
        print $out_fh $_."\n";
        next;
    }
    
    my @a = split /\s+/;
    my @alt_alleles = split /,/, $a[4];
    my $ref_len = length($a[3]);
    my $num_SV = 0;
    
    for my $alt (@alt_alleles) {
        my $alt_len = length($alt);
        if (abs($alt_len - $ref_len) >= 50) {
            $num_SV++;
        }
    }
    
    $total++;
    
    if ($num_SV > 0) {
        print $out_fh $_."\n";
        $sv_count++;
    }
}

close $in_fh;
close $out_fh;

print STDERR "\nSummary:\n";
print STDERR "  Total variants: $total\n";
print STDERR "  SV variants:    $sv_count\n";
print STDERR "  Non-SV:         ", $total - $sv_count, "\n";
print STDERR "  Output:         $out_file\n";
