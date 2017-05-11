#!/usr/bin/env perl

use strict;
use warnings;

my $wantedmiF = $ARGV[0];
my $outF = $ARGV[1];

my $line;
#-------------------------------
open (IN,$wantedmiF) or die;
open OUT,">".$outF or die;
open OUTCS,">".$outF."_cs" or die;

while ($line=<IN>)
{
  chomp $line;
  if ($line =~ m/^>/)
  {
    $line =~ s/>//;
    my $mi = $line;
    $line = <IN>;	chomp $line;
    	my $seq = uc($line);

  print OUT $mi,"\t",substr($seq,1,7),"\t9606\n"; 
  print OUTCS $mi,"\t9606\t",$mi,"\t",$seq,"\n";
  }
}

close IN;
close OUT;	close OUTCS;


