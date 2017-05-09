#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in-dir> <trackname> <out-dir>\n" unless @ARGV==3;
my ($indir,$trackName,$outdir)=@ARGV;

my @files=`ls $indir`;
foreach my $file (@files) {
  chomp $file;
  my $found=0;
  my $listening=0;
  open(IN,"$indir/$file") || die $file;
  while(<IN>) {
    if(/^%(\S+)/) {$listening=$1 eq $trackName; next}
    elsif(/[1-9]/){if($listening) {$found=1; last } }
  }
  close(IN);
  if($found) { system("cp $indir/$file $outdir") }
}
