#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $NUM_BINS=300;

my $name=ProgramName::get();
die "$name <track-name> <fastb-dir> | histogram\n" unless @ARGV==2;
my ($trackName,$dir)=@ARGV;

my (%counts,$N);
my @files=`ls $dir`;
foreach my $file (@files) {
  my $echoing=0;
  open(IN,"$dir/$file") || die "can't open file: $dir/$file\n";
  while(<IN>) {
    chomp;
    if(/>/) {$echoing=0}
    elsif(/%(\S+)/) {$echoing=($1 eq $trackName)}
    elsif($echoing && $_>0) {
      my $bin=int($_*$NUM_BINS);
      ++$counts{$bin};
      ++$N;
      #print
    }
  }
  close(IN);
}
my $cum=0;
my $firstVal=1/6;
my $secondVal=3/6;
my $thirdVal=5/6;
my (%sums,%sampleSizes);
for(my $bin=0 ; $bin<=$NUM_BINS ; ++$bin) {
  my $P=$counts{$bin}/$N;
  my $newCum=$cum+$P;
  my $segment;
  if($newCum<$firstVal) {$segment=1}
  elsif($newCum<$secondVal) {$segment=2}
  else {$segment=3};
  $sums{$segment}+=$P*$bin/$NUM_BINS;
  $sampleSizes{$segment}+=$P;
  #if($cum<$firstVal && $newCum>=$firstVal ||
  #   $cum<$secondVal && $newCum>=$secondVal ||
  #   $cum<$thirdVal && $newCum>=$thirdVal) {
  #  my $val=$bin/$NUM_BINS;
  #  print "$newCum\t$val\n";
  #}
  #print "$bin\t$P\n";
  $cum+=$P
}
for(my $i=1 ; $i<=3 ; ++$i) {
  my $y=$sums{$i}/$sampleSizes{$i};
  print "$i\t$y\n";
}
