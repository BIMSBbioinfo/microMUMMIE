#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <*.fastb>\n" unless @ARGV==1;
my ($infile)=@ARGV;

open(IN,$infile) || die;
while(<IN>) {
  if(/^\%\s*(\S+)/) {print "$1 : continuous\n"}
  elsif(/^>\s*(\S+)/) {print "$1 : discrete ACGT\n"}
}
close(IN);

