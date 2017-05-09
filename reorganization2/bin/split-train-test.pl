#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <indir> <num-partitions> <dir-stem>\n" unless @ARGV==3;
my ($indir,$numPartitions,$dirStem)=@ARGV;

for(my $i=1 ; $i<=$numPartitions ; ++$i) {
  if(-e "$dirStem-train-$i") { die "$dirStem-train-$i already exists\n" }
  if(-e "$dirStem-test-$i") { die "$dirStem-test-$i already exists\n" }
  system("mkdir $dirStem-train-$i");
  system("mkdir $dirStem-test-$i");
}

my @infiles=`ls $indir | grep .fastb | shuffle`;
my $n=@infiles;
my $partitionSize=int($n/$numPartitions);
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=@infiles[$i];
  chomp $file;
  my $partition=int($i/$partitionSize)+1;
  if($partition>$numPartitions) {$partition=$numPartitions}
  my $dir="$dirStem-test-$partition";
  system("ln -s ../$indir/$file $dir/$file");
  for(my $j=1 ; $j<=$numPartitions ; ++$j) {
    if($j!=$partition) {
      my $dir="$dirStem-train-$j";
      system("ln -s ../$indir/$file $dir/$file");
    }
  }
}


