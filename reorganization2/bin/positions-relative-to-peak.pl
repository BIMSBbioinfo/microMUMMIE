#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <peaks.gff> <sites.gff>\n" unless @ARGV==2;
my ($peaksFile,$sitesFile)=@ARGV;

my %peaks;
open(IN,$peaksFile) || die;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my ($substrate,$x,$y,$peak)=@fields;
  push @{$peaks{$substrate}},$peak;
}
close(IN);

open(IN,$sitesFile) || die;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my ($substrate,$x,$y,$begin,$end)=@fields;
  my $peaks=$peaks{$substrate};
  my ($closestPeak,$closestDistance);
  foreach my $peak (@$peaks) {
    my $dist=min(abs($peak-$begin),abs($peak-$end));
    if(!defined($closestPeak) || $dist<$closestDistance) {
      $closestPeak=$peak;
      $closestDistance=$dist;
    }
  }
  my $pos=$begin-$closestPeak;
  print "$pos\n";
}
close(IN);

sub min {
  my ($a,$b)=@_;
  return $a<$b ? $a : $b;
}


