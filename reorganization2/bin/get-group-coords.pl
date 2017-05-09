#!/usr/bin/perl
use strict;

die "get-group-coords.pl <in.groups> <out.txt>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $nextID=1;
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
open(IN,$infile) || die "can't open file: $infile\n";
<IN>;
while(<IN>) {
  chomp;
  my @fields=split/,/,$_;
  my ($chr,$strand,$begin,$end)=@fields;
  print OUT "$chr\t$begin\t$end\tGENE$nextID\t$strand\tENST0000$nextID\n";

  ++$nextID;
}
close(IN);
close(OUT);





