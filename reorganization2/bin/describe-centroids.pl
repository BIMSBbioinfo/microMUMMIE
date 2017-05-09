#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <centroids.txt> <schema.txt>\n" unless @ARGV==2;
my ($centroidsFile,$schemaFile)=@ARGV;

my @features;
open(IN,$schemaFile) || die "can't open $schemaFile\n";
while(<IN>) {
  if(/^([^\s:=]+)/) {push @features,$1}
}
close(IN);

open(IN,$centroidsFile) || die "can't open $centroidsFile\n";
my $N=0+<IN>;
my $comp=1;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my $n=0+shift @fields;
  my %hash;
  for(my $i=0 ; $i<$n ; ++$i) { $hash{$features[$i]}=$fields[$i]  }
  print "\ncomponent #$comp: ";
  my @sorted=sort {$hash{$b} <=> $hash{$a}} @features;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $y=$sorted[$i];
    print "$y\t";
  }
  print "\n";
  ++$comp;
}
close(IN);




