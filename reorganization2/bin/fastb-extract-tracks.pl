#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.fastb> <out.fastb> <track1> <track2> ...\n" unless @ARGV>=3;
my $infile=shift @ARGV;
my $outfile=shift @ARGV;
my %keep;
while(@ARGV) {
  my $track=shift @ARGV;
  $keep{$track}=1;
}

open(OUT,">$outfile") || die "can't create file $outfile\n";
open(IN,$infile) || die "can't open file $infile\n";
my $echo=0;
while(<IN>) {
  if(/^[%>](\S+)/) { $echo=$keep{$1} }
  else {$_=~s/nan/0/g}
  if($echo) { print OUT }
}
close(IN);
close(OUT);

