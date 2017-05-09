#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <track1> <track2> <in-dir1> <in-dir2> <out-dir>\n"
  unless @ARGV==5;
my ($trackName1,$trackName2,$fastbDir1,$fastbDir2,$outDir)=@ARGV;

my @files=`ls $fastbDir1`;
my $numFiles=@files;
my $errors;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  $file=~/([^\/]+)\.fastb/ || die;
  my $stem=$1;
  my $fastb1=new Fastb("$fastbDir1/$file");
  if(-e "$fastbDir2/$file") {
    my $fastb2=new Fastb("$fastbDir2/$file");
    my $track1=$fastb1->getTrackByName($trackName1);
    my $track2=$fastb2->getTrackByName($trackName2);
    die unless $track1 && $track2;
    my $data1=$track1->getData();
    my $data2=$track2->getData();
    my $L1=@$data1;
    my $L2=@$data2;
    #die "unequal lengths in file $file: $L1 vs. $L2\n" unless $L1==$L2;
    if($L1!=$L2) { ++$errors; next }
    my $intervals1=$track1->getNonzeroRegions();
    my $intervals2=$track2->getNonzeroRegions();
    foreach my $interval1 (@$intervals1) {
      foreach my $interval2 (@$intervals2) {
	if(overlap($interval2,$interval1)) {
	  my ($begin,$end)=@$interval1;
	  for(my $i=$begin ; $i<$end ; ++$i) { $data1->[$i]=0 }
	}
      }
    }
    my $nonzero=0;
    for(my $i=0 ; $i<$L1 ; ++$i)
      { if($data1->[$i]>0) { $nonzero=1; last } }
    next unless $nonzero;
  }
  $fastb1->save("$outDir/$file");
}
print "$errors errors / $numFiles files\n";

sub overlap {
  my ($a,$b)=@_;
  return $a->[0]<$b->[1] && $b->[0]<$a->[1];
}


