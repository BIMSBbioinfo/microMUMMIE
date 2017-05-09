#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my %hash;
my @files=`ls trained*.hmm`;
foreach my $file (@files) {
  chomp $file;
  my $msg=`MUMMIE/summarize-hmm $file`;
  $msg=~/state (\d+) is probably foreground/ || die $msg;
  my $fg=$1;
  $msg=~/state (\d+) is probably background/ || die $msg;
  my $bg=$1;
  my @msg=`MUMMIE/extract-motifs $file 10 8 -s $fg -L $bg`;
  foreach my $msg (@msg) {
    $msg=~/^(\S+)/;
    ++$hash{$1};
  }
}
my @keys=keys %hash;
@keys=sort {$hash{$b} <=> $hash{$a}} @keys;
foreach my $key (@keys) {
  my $count=$hash{$key};
  print "$key\t$count\n";
}

