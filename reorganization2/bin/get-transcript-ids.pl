#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <fastb-dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

my @files=`ls $dir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i];
  next unless $file=~/([^\/]+).fastb/;
  my $name=$1;
  open(IN,"$dir/$file") || die "can't open $dir/$file\n";
  while(<IN>) {
    if(/transcriptID=(\S+)/) { print "$name\t$1\n"; last }
  }
  close(IN);
}


