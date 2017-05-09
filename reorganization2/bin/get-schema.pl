#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

die "get-schema.pl <fastb-dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;
my @files=`ls $dir`;
foreach my $file (@files) {
  chomp $file;
  if($file=~/\.fastb$/) {
    system("schema-from-fastb.pl $dir/$file");
    exit;
  }
}
