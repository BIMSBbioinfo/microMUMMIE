#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use Fastb;

my $WANT_MASKING=1;

my $name=ProgramName::get();
die "$name <track1> <track2> <in-dir1> <in-dir2> <out-dir> <threshold>\n" 
  unless @ARGV==6;
my ($trackName1,$trackName2,$fastbDir1,$fastbDir2,$outDir,$threshold)=@ARGV;

my @files=`ls $fastbDir1`;
my $numFiles=@files;
my $errors=0;
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
    my $anyNonzero=0;
    for(my $i=0 ; $i<$L1 ; ++$i) {
      my $a=$data1->[$i];
      my $b=$data2->[$i];
      if($WANT_MASKING && $b>0) { $b=$a }
      $a-=$b;
      if($a<0) { $a=0 }
      $data1->[$i]=$a;
      #if($a>$threshold) { $nonzero=1 }
    }
    my $intervals=$track1->getNonzeroRegions();
    foreach my $interval (@$intervals) {
      my $nonzero=0;
      my ($begin,$end)=@$interval;
      for(my $i=$begin ; $i<$end ; ++$i)
	{ if($data1->[$i]>$threshold) { $nonzero=1; last } }
      if($end-$begin<6) { $nonzero=0 }
      if($nonzero) { $anyNonzero=1 }
      else { for(my $i=$begin ; $i<$end ; ++$i) { $data1->[$i]=0 } }
    }
    next unless $anyNonzero;
  }
  $fastb1->save("$outDir/$file");
}
print "$errors errors / $numFiles files\n";


