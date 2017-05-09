#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <schema.txt> <full|diag> <outfile.tgf>\n" unless @ARGV==3;
my ($schemaFile,$fullDiag,$outfile)=@ARGV;
$fullDiag="\U$fullDiag";
my $full;
if($fullDiag eq "FULL") {$full=1}
elsif($fullDiag eq "DIAG") {$full=0}
else {die "specify FULL or DIAG for covariance matrix\n"}

open(OUT,">$outfile") || die "can't write to file $outfile\n";
my $id=1;
open(IN,$schemaFile) || die "Can't open $schemaFile\n";
while(<IN>) {
  if(/(\S+)\s*:\s*continuous/) {
    print OUT "$id $1\n";
    ++$id;
  }
}
close(IN);
print OUT "#\n";
my $n=$id-1;
for(my $i=1 ; $i<=$n ; ++$i) {
  for(my $j=1 ; $j<=$n ; ++$j) {
    print OUT "$i\t$j\n";
  }
}
close(OUT);



