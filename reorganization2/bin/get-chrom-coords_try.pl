#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff> \n" unless @ARGV==1;
my ($infile)=@ARGV;

open(IN,$infile) || die "can't open file: $infile\n";
while(<IN>) {
  chomp;
  push my @data, $_;
  print @data,"\n";
}
close(IN);

