#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.hmm> <out.hmm>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die;
open(IN,$infile) || die;
while(<IN>) {
  if(/transitions/) { print OUT; last }
  else {print OUT}
}
$_=<IN>; # table size
print OUT;
$_=<IN>; my @fields=split/\s+/,$_; # state 0
my $epsilon=$fields[1];
print OUT;
$_=<IN>; @fields=split/\s+/,$_; # state 1 (background)
$epsilon+=$fields[1];
print OUT;
$_=<IN>; my @fieldsA=split/\s+/,$_; # state 2 (factor A)
my $A=$fieldsA[1];
$_=<IN>; my @fieldsB=split/\s+/,$_; # state 3 (factor B)
my $B=$fieldsB[1];
$_=<IN>; my @fieldsC=split/\s+/,$_; # state 4 (both factors)
my $r=$A/$B;
my $newB=(-1-$r+sqrt(($r+1)*($r+1)-(4*$r*($epsilon-1))))/(2*$r);
my $newA=$r*$newB;
my $newC=$newA*$newB;
$fieldsA[1]=$newA;
$fieldsB[1]=$newB;
$fieldsC[1]=$newC;
my $line=join("\t",@fieldsA); print OUT "$line\n";
my $line=join("\t",@fieldsB); print OUT "$line\n";
my $line=join("\t",@fieldsC); print OUT "$line\n";
while(<IN>) {print OUT}
close(IN);
close(OUT);
