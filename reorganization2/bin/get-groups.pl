#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use Fastb;
use FastbTrack;

my $name=ProgramName::get();
die "$name <main-track> <fastb-dir> <extra-margin>\n" unless @ARGV==3;
my ($mainTrackName,$fastbDir,$margin)=@ARGV;

my @files=`ls $fastbDir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  $file=~/([^\/]+)\.fastb/ || next;
  my $stem=$1;
  my $fastb=new Fastb("$fastbDir/$file");
  my $track=$fastb->getTrackByName($mainTrackName);
  die "can't find track: $mainTrackName\n" unless $track;
  my $L=$track->getLength();
print "Length=$L","\n";
  my $array=$track->getNonzeroRegions();
  foreach my $interval (@$array) {
    my ($begin,$end)=@$interval;
    ++$begin;
    my $b=$begin-$margin;
    my $e=$end+$margin;
    if($b<1) {$b=1}
    if($e>$L) {$e=$L}
    print "$stem\tPAR\tgroup\t$b\t$e\t.\t+\t.\n";
  }
}



