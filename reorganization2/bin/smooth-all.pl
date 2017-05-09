#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in-dir> <out-dir> <window-size> <#iterations> [tracks]\n" unless @ARGV>=4;
my ($inDir,$outDir,$window,$iter,$trackList)=@ARGV;

if(!-e $outDir) {system("mkdir $outDir")}
my $dashT=length($trackList)>0 ? "-t $trackList" : "";

my @files=`ls $inDir`;
foreach my $file (@files) {
  chomp $file;
  if($file=~/([^\/]+)\.fastb$/) {
    my $stem=$1;
    my $outfile="$outDir/$stem.fastb";
    my $cmd="smooth-fastb $dashT $inDir/$file $window $iter $outfile";
    system($cmd);
  }
}




