#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <old-track-name> <new-track-name> <fastb-dir>\n" unless @ARGV==3;
my ($oldName,$newName,$dir)=@ARGV;

my @files=`ls $dir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  next unless $file=~/\.fastb/;
  my $path="$dir/$file";
  my $fastb=new Fastb($path);
  my $track=$fastb->renameTrack($oldName,$newName);
  $fastb->save($path);
}



