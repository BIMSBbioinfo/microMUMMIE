#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.fastb> <out.gff> <trackname> <threshold>\n" unless @ARGV==4;
my ($infile,$outfile,$wantTrack,$threshold)=@ARGV;

$infile=~/([^\/]+)\.fastb/ || die;
my $substrate=$1;
open(OUT,">$outfile") || die;
open(IN,$infile) || die;
my $x=0;
my $echo=0;
my $inSite=0;
my ($begin,$end);
while(<IN>) {
  if(/[%>]\s*(\S+)/) {
    if($1 eq $wantTrack) { $echo=1 } else { $echo=0 }
    $x=0;
  }
  elsif(/(\S+)/) {
    if($echo) {
      my $y=0+$_;
      if($y>=$threshold) {
	if(!$inSite) { $inSite=1; $begin=$x+1; }
      }
      else {
	if($inSite) {
	  $inSite=0;
	  $end=$x+1;
	  print OUT "$substrate\t$wantTrack\tsite\t$begin\t$end\t.\t+\t.\n";
	}
      }
    }
    ++$x;
  }
}
close(IN);
close(OUT);
