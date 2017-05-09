#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <dir> <num-bins> <states>\n" unless @ARGV==3;
my ($dir,$numBins,$keepStates)=@ARGV;

my %keepStates;
my @fields=split/[,\s]/,$keepStates;
foreach my $q (@fields) {$keepStates{$q}=1}

my (@counts,%states);
my @files=`ls $dir`;
foreach my $file (@files) {
  chomp $file;
  open(IN,"$dir/$file") || die "can't open $dir/$file\n";
  my @path;
  while(<IN>) {
    if(/(\d+)/){push @path,$1}
  }
  close(IN);
  my $L=@path;
  my $binSize=$L/$numBins;
  for(my $i=0 ; $i<$L ; ++$i) {
    my $q=$path[$i];
    ++$states{$q};
    my $bin=int($i/$binSize);
    if($bin>=$numBins) {$bin=$numBins-1}
    ++$counts[$bin]->[$q];
  }
}

my @states=keys %states;
@states=sort {$a <=> $b} @states;
my @binTotals;
for(my $i=0 ; $i<$numBins ; ++$i) {
  my $bin=$counts[$i];
  foreach my $q (@states) {
    $binTotals[$i]+=$bin->[$q];
  }	
}
foreach my $q (@states) {
  next unless($keepStates{$q});
  for(my $i=0 ; $i<$numBins ; ++$i) {
    my $P=$counts[$i]->[$q]/$binTotals[$i];
    print "$i\t$P\n";
  }
  print "\n";
}

