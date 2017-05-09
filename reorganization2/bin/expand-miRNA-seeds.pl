#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use FastaReader;
use Fastb;

my $name=ProgramName::get();
die "$name <8mers.fasta>\n" unless @ARGV==1;
my ($motifFile)=@ARGV;

my $reader=new FastaReader($motifFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/>(\S+)/ || die;
  print substr($sequence,0,8) . "\t$1" . "\n";      # 8mer-m1
  print substr($sequence,0,7)."A" . "\t$1" . "\n";  # 8mer-A1
  print substr($sequence,1,7) . "\t$1" . "\n";      # 7mer-m1
  print substr($sequence,1,6)."A" . "\t$1" . "\n";  # 7mer-A1
  print substr($sequence,0,7) . "\t$1" . "\n";      # 7mer-m8
  print substr($sequence,0,6) . "\t$1" . "\n";      # 6mer3-8
  print substr($sequence,1,6) . "\t$1" . "\n";      # 6mer1-7
}
