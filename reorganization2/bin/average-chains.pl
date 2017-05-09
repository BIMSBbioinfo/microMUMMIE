#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <chain1> <chain2> <chain3> ...\n" unless @ARGV>1;
my (@chains,$ORDER,$ALPHA);
foreach my $filename (@ARGV) {
  push @chains,load($filename);
}
my $numChains=@chains;
print "$ORDER order\nalphabet:\n$ALPHA\n";
#my $ave=[];
my $L=@{$chains[0]};
for(my $j=0 ; $j<$L ; ++$j) {
  my $sum;
  my $nmer;
  for(my $i=0 ; $i<$numChains ; ++$i) {
    my $chain=$chains[$i];
    my $elem=$chain->[$j];
    if(!$nmer) {$nmer=$elem->[0]}
    if($nmer && $nmer ne $elem->[0]) { die "$nmer != $elem->[0]"}
    $sum+=exp($elem->[1])/$numChains;
  }
  #push @$ave,[$nmer,$sum];
  my $val=log($sum);
  print "$nmer\t$val\n";
}

sub load
  {
    my ($filename)=@_;
    my $chain=[];
    open(IN,$filename) || die "Can't open file: $filename\n";
    while(<IN>) {
      if(/(\d+)\s+order/) {$ORDER=$1}
      elsif(/alphabet:/) {
	$ALPHA=<IN>;
	chomp $ALPHA;
      }
      elsif(/(\S+)\s+(\S+)/) { push @$chain,[$1,$2] }
    }
    close(IN);
    return $chain;
  }




