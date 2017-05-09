#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $CONTEXT_MARGIN=200;
my $INTERVAL_SIZE=8; # the -u option to RNAplfold

#============================================================
$CONTEXT_MARGIN-=($INTERVAL_SIZE-1);
my $name=ProgramName::get();
die "$name <indir-fastb> <indir-fold> <outdir>\n" unless @ARGV==3;
my ($fastbDir,$foldDir,$outdir)=@ARGV;
my @fastbFiles=`ls $fastbDir`;
foreach my $fastb (@fastbFiles) {
  $fastb=~/([^\/]+)\.fastb/;
  my $stem=$1;
  my $fastbFile="$fastbDir/$stem.fastb";
  open(FASTB,$fastbFile) || die "can't open $fastbFile\n";
  my $outfile="$outdir/$stem.fastb";
  open(OUT,">$outfile") || die "can't write to file $outfile\n";
  my ($transcriptID,$fastbLength);
  while(<FASTB>) {
    if(/transcriptID=(\S+)/) {$transcriptID=$1}
    if(/length=(\d+)/) {$fastbLength=$1}
    print OUT;
  }
  close(FASTB);
  #my $foldFile="$foldDir/$stem.fold";
  my $foldFile="$foldDir/$transcriptID.fold";
  if(!open(FOLD,$foldFile)) {# || die "can't open $foldFile\n";
    close(OUT);
    system("rm $outfile");
    next;
  }
  for(my $i=0 ; $i<2 ; ++$i) { <FOLD> }
  #for(my $i=0 ; $i<$CONTEXT_MARGIN ; ++$i) { <FOLD> }
  my @F;
  while(<FOLD>) {
    if(/(\S+)\s+(\S+)/) { push @F,$2 }
  }
  while(@F>$fastbLength) {shift @F}
  print OUT "%unpaired\n";
  my $FN=@F;
  for(my $i=0 ; $i<$FN ; ++$i) {
    my $x=$F[$i];
    print OUT "$x\n";
  }
  close(OUT);
  close(FOLD);
}




