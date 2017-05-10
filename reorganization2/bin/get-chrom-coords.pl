#!/usr/bin/env perl
use strict;
use ProgramName;



my $name=ProgramName::get();
die "$name <in.gff> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(IN,$infile) || die "can't open file: $infile\n";
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
while(<IN>) {
  chomp;
  my @fields=split/\t/,$_;
  my $substrate=$fields[0];
  my $begin=$fields[3];
  my $end=$fields[4];
  my $extra=$fields[8];
 my ($seq,$type,$miRNA);
  if($extra=~/seq=([^;]+);\s*type=([^;]+);\s*miRNA=([^;]+);/)
    { ($seq,$type,$miRNA)=($1,$2,$3); }
  else { $extra=~/seq=([^;]+);\s*miRNA=([^;]+);/ || die;
	 ($seq,$type,$miRNA)=($1,"site",$2); }
  $fields[1]=$miRNA;
  $fields[2]=$type;
  $fields[8]="seq=$seq;";
  $substrate=~ /(\w+).fastb/ ||
    die "can't parse substrate: $substrate\n";
  my ($chr,$chunkBegin,$strand,$chunkEnd)=($1,$2,$3,$4);
  $fields[6]=$strand;
  if($strand eq "+") {
    $begin+=$chunkBegin;
    $end+=$chunkBegin;
  }
  else {
    my $b=$begin; my $e=$end; my $len=$e-$b;
    $begin=$chunkEnd-$end+1;
    $end=$begin+$len;
  }
  $fields[0]=$chr;
  $fields[3]=$begin;
  $fields[4]=$end;
  my $line=join("\t",@fields);
  print OUT "$line\n";
}
close(OUT);
close(IN);

