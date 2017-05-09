#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use FastaReader;
use FastaWriter;

my $name=ProgramName::get();
die "$name <in.fasta> <in.schema> <out-dir>\n" unless @ARGV==3;
my ($fasta,$schemaFile,$outDir)=@ARGV;

my (@continuous,@discrete);
open(IN,$schemaFile) || die "can't open $schemaFile\n";
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my $trackName=$fields[0];
  my $type=$fields[2];
  my $alphabet=$fields[3];
  $type="\L$type";
  if($type eq "continuous") { push @continuous,$trackName }
  else { push @discrete,$trackName }
}
close(IN);
if(@discrete>1) {die "schema contains multiple discrete tracks"}
my $discreteTrack=$discrete[0];

my $reader=new FastaReader($fasta);
my $writer=new FastaWriter;
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/^\s*>\s*(\S+)/ || die "can't parse defline: $defline\n";
  my $id=$1;
  my $filename="$outDir/$id.fastb";
  open(OUT,">$filename") || die "can't write to file: $filename\n";
  #print OUT "$defline\n$sequence\n";
  $writer->addToFasta(">$discreteTrack",$sequence,\*OUT);
  my $L=length($sequence);
  foreach my $cont (@continuous) {
    print OUT "%$cont\n";
    for(my $i=0 ; $i<$L ; ++$i) { print OUT "0\n" }
  }
  close(OUT);
}






