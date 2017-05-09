#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffReader;
my $features=$reader->loadGFF($infile);
my $n=@$features;
for(my $i=0 ; $i+1<$n ; ++$i) {
  my $this=$features->[$i];
  my $next=$features->[$i+1];
  if($this->getBegin()==$next->getBegin() &&
     $this->getEnd()==$next->getEnd()) {
    my $thisLen=$this->getLength();
    my $nextLen=$next->getLength();
    my $thisScore=$this->getScore();
    my $nextScore=$next->getScore();
    my $combinedScore=$thisScore+$nextScore;
    my $best=$thisLen>=$nextLen ? $this : $next;
    $best->setScore($combinedScore);
    splice(@$features,$i+1,1);
    $features->[$i]=$best;
    --$n;
    --$i;
  }
}
for(my $i=0 ; $i<$n ; ++$i) {
  my $feature=$features->[$i];
  print $feature->toGff();
}




