#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.fastb> <out.fasta>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

$infile=~/([^\/]+)\.fastb/ || die "can't parse filename\n";
my $stem=$1;

open(OUT,">$outfile") || die "can't create file $outfile\n";
open(IN,$infile) || die "can't open file $infile\n";
my $echo=0;
while(<IN>) {
  if(/^%/) { $echo=0 }
  elsif(/^>\s*(\S+)/) { $echo=1; print OUT ">$stem /track=$1\n"}
  elsif($echo) { print OUT }
}
close(IN);
close(OUT);

