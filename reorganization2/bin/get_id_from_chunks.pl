#!/bin/perl
use strict; use warnings;

my $defidF="AGO_defid.txt";
my $outF = "AGO_enst.txt";

open IN,$defidF or die;
open OUT,">".$outF or die;

while (<IN>)
{
  if ($_ =~ m/(ENST\d+) /) { print OUT $1,"\n"; }
}

close IN; close OUT;

