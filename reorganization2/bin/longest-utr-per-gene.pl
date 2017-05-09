#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <in-dir> <out-dir>\n" unless @ARGV==2;
my ($inDir,$outDir)=@ARGV;

my @files=`ls $inDir`;
@files=sort {$a<=>$b} @files;
my %lengths;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  my $fastb=new Fastb("$inDir/$file");
  my $L=$fastb->getLength();
  $file=~/([^_]+)/ || die $file;
  my $gene=$1;
  $lengths{$gene}->{$file}=$L;
}

my @geneNames=keys %lengths;
my $numGenes=@geneNames;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $geneName=$geneNames[$i];
  my $gene=$lengths{$geneName};
  my @files=keys %$gene;
  @files=sort {$gene->{$b}<=>$gene->{$a}} @files;
  my $longest=$files[0];
  chomp $longest;
  system("cp $inDir/$longest $outDir");
}






