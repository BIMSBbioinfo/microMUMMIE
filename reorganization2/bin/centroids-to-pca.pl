#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

open(DATA,">pca.data");
my $first=1;
<STDIN>;
while(<STDIN>) {
  chomp;
  next unless($_=~/\S/);
  my @fields=split/\s+/,$_;
  shift @fields;
  my $line=join(" ",@fields);
  unless($line=~/inf/ || $line=~/nan/) { print DATA "$line 0\n" }
  if($first) {
    open(NAMES,">pca.names");
    print NAMES "2 categories\n";
    my $n=@fields;
    for(my $x=0 ; $x<$n ; ++$x) {
      print NAMES "x$x: continuous\n";
    }
    close(NAMES);
  }
}
close(DATA);
