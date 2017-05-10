#!/usr/bin/env perl
use strict;

my %peaks;
open(IN,"peaks.gff") || die;
while(<IN>) {
    chomp;
    my @fields=split/\t/,$_;
    my ($chunkID,$program,$elemType,$begin,$end,$score,$strand)=@fields;
    push @{$peaks{$chunkID}}, $begin;
}
close(IN);

open(IN,"output.gff") || die;
while(<IN>) {
    chomp;
    my @fields=split/\t/,$_;
    my ($chunkID,$program,$elemType,$begin,$end,$score,$strand)=@fields;
    my $peaks=$peaks{$chunkID};
    my ($closest,$bestDist);
    foreach my $peak (@$peaks) {
	my $dist;
	if($begin>=$peak) { $dist=$begin-$peak }
	elsif($end<=$peak) { $dist=$peak-$end }
	else { $dist=0 }
	if($dist<$bestDist || !defined($bestDist)) {
	    $bestDist=$dist;
	    $closest=$peak;
	}
    }
    my $relPos=$begin-$closest;
    print "$relPos\n";
}
close(IN);


