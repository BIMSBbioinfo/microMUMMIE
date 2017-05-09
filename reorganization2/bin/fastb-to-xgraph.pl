#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my $LINEWIDTH=3;

my $wantTrack;
if(@ARGV>0) {$wantTrack=shift @ARGV}

my (@names,$stdout);
my $x=0;
my $trackID=0;
my $echo=0;
while(<STDIN>) {
  if(/%\s*(\S+)/) {
    ++$trackID;
    if(!defined($wantTrack) || $1 eq $wantTrack) {
      push @names,$1;
      $stdout.="\n";
      $x=0;
      $echo=1;
    }
    else { $echo=0 }
  }
  elsif(/^\s*>/) {$echo=0}
  elsif(/(\S+)/) {
    next unless $echo;
    $stdout.="$x $1\n";
    ++$x;
  }
}
my $names="";
my $n=@names;
for(my $i=0 ; $i<$n ; ++$i) {$names.="-$i $names[$i] "}
open(OUT,"|xgraph -lw $LINEWIDTH $names");
print OUT $stdout;
close(OUT);
#print $stdout;
