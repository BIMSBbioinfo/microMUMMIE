#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <clusters.csv>\n" unless @ARGV==1;
my ($infile)=@ARGV;

open(IN,$infile) || die "can't open $infile\n";
while(<IN>) {
  my @fields=split/,/,$_;
  my ($Chromosome,$Strand,$ClusterStart,$ClusterEnd,$ClusterID,$ClusterSequence,$ReadCount,$ModeLocation,$ModeScore,$ConversionLocationCount,$ConversionEventCount,$NonConversionEventCount,$FilterType,$TranscriptLocation,$TranscriptID,$GeneName)=@fields;
  next unless $TranscriptLocation eq "3UTR";
  my $site=
    {
     chr=>$Chromosome,
     strand=>$Strand,
     begin=>$ClusterStart,
     end=>$ClusterEnd,
     gene=>$GeneName,
     transcript=>$TranscriptID
    };
}
close(IN);

