#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use Translation;
use NgramIterator;
use FastaReader;

my $SHOULD_LOGIFY_WEIGHTS=1;
my $CONTEXT_LENGTH=0;

my $name=ProgramName::get();
die "$name <mature-mirnas.fasta> <threshold> <weights.out> <sequences.out>\n"
  unless @ARGV==4;
my ($infile,$threshold,$weightsFile,$outfile)=@ARGV;

my %seedsSeen;
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
open(WEIGHT,">$weightsFile") || die "can't write to file: $weightsFile\n";
my $reader=new FastaReader($infile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  my $weight=1;
  $defline=~/>(\S+)/;
  my $ID=$1;
  next if($ID=~/hsa-mir-3187/ || $ID=~/hsa-mir-636/);

  if($defline=~/\/weight=(\S+)/) {$weight=$1}
  next unless $weight>$threshold;
  $sequence=~s/U/T/g;
  my $eightmer=substr($sequence,0,8);
  my $revSeq=Translation::reverseComplement(\$eightmer);
  my $iter=new NgramIterator("ATCG",$CONTEXT_LENGTH);
  if($CONTEXT_LENGTH>0) {
    while(1) {
      my $context=$iter->nextString();
      last unless $context;
      emit("$context$revSeq",$ID,$weight);
    }
  }
  else {
    emit($revSeq,$ID,$weight);
  }
}
$reader->close();
close(OUT);
close(WEIGHT);

sub emit
  {
    my ($string,$ID,$weight)=@_;
    if(!$seedsSeen{$string}) {
      print OUT ">$ID\n$string\n";
      if($SHOULD_LOGIFY_WEIGHTS) { $weight=log($weight) }
      print WEIGHT "$ID\t$weight\n";
      $seedsSeen{$string}=1;
    }
  }
