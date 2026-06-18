#!/usr/bin/perl -w
use strict;
use warnings;
my $file=shift;
open IN, "gzip -cd  $file |" or die $!;

while (<IN>) {
    chomp;
    if (/^#/){
    print $_."\n";
    next;
    }
    my @a = split/\s+/;
    my @b = split(/,/, $a[4]);
    my $k = 0;
    my @d = split(/;/, $a[7]);
    $d[0] =~ s/AC=//g;
    my @e = split(/,/, $d[0]);
    my $len1 = length($a[3]);
    my $num_SV=0;
    my $num_SNP;
    my @filtered_e;
    my $sum=0;
    for (my $i = 0; $i < @b; $i++) {
        my $len_q = length($b[$i]);
        if (abs($len_q -$len1 )>= 50) {
            $num_SV++;
        }
    }
    if ($num_SV > 0) {
       print  $_."\n";
    }
}

close IN;
