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
die "$name <mature-mirnas.txt> <outfile.fasta>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my %hash;
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
open(IN,$infile) || die "can't open $infile\n";
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my ($ID,$sequence)=@fields;
  #next if($ID=~/hsa-mir-3187/ || $ID=~/hsa-mir-636/);
  $sequence=~s/U/T/g;
  my $eightmer=substr($sequence,0,8);
  my $revSeq=Translation::reverseComplement(\$eightmer);
  my $iter=new NgramIterator("ATCG",$CONTEXT_LENGTH);
  #emit($revSeq,$ID);
  push @{$hash{$revSeq}},$ID;
}
close(IN);

my @keys=keys %hash;
my $numKeys=@keys;
for(my $i=0 ; $i<$numKeys ; ++$i) {
  my $key=$keys[$i];
  my $IDs=$hash{$key};
  my $numIDs=@$IDs;
  my $index=int(rand($numIDs));
  emit($key,$IDs->[$index]);
}
close(OUT);


my %seedsSeen;
sub emit
  {
    my ($string,$ID)=@_;
    if(!$seedsSeen{$string}) {
      print OUT ">$ID\n$string\n";
      $seedsSeen{$string}=1;
    }
  }

