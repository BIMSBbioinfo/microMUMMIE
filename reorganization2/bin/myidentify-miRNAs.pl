#!/usr/bin/env perl
use strict;
use GffReader;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <in.gff> <seeds.fasta>\n" unless @ARGV==2;
my ($gffFile,$seedsFile)=@ARGV;

my (%motifs);
my $reader=new FastaReader($seedsFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/>(\S+)/ || die;

   $motifs{$sequence}->{$1}=1;      # exact bulge sm
   print $sequence,"\t",$1,"\n";
}

my $reader=new GffReader();
my $sites=$reader->loadGFF($gffFile);

my $n=@$sites;
for(my $i=0 ; $i<$n ; ++$i) {
  my $site=$sites->[$i];
  my $seq;
  my $extra=$site->{additionalFields};
  my @path;
  my $numExtra=@$extra;
  my $seqField=-1;
  for(my $i=0 ; $i<$numExtra ; ++$i) {
    my $field=$extra->[$i];
    if($field=~/seq=([^;]+);/) { $seq=$1; $seqField=$i }
  print $seq,"\n",$seqField,"\n"; die;
    #if($field=~/path=([^;]+)/) { @path=split/,/,$1 }
  }
  #die unless @path>0;
  my %hits;
  addHits(\%hits,\%motifs,$seq);

 # if(@path==0) {
 #   if($INCLUDE_8MERS)  { addHits(\%hits,\%motifs1_8,substr($seq,0,8)) }
  #  if($INCLUDE_7MERS)  { addHits(\%hits,\%motifs1_7,substr($seq,0,7)) }
 #   if($INCLUDE_7MERS)  { addHits(\%hits,\%motifs2_8,substr($seq,1,7)) }
 #   if($INCLUDE_6MERS)  { addHits(\%hits,\%motifs2_7,substr($seq,1,6)) }
 #    if($INCLUDE_6MERS)  { addHits(\%hits,\%motifs1_6,substr($seq,0,6)) }
 # }

  my @miRNAs=keys %hits;
  my @extra=@$extra;
  foreach my $miRNA (@miRNAs) {
    my @ex=@extra;
    push @ex,"miRNA=$miRNA;";
    $site->{additionalFields}=\@ex;
    my $gff=$site->toGff();
    print "$gff";
  }
}


sub addHits {
  my ($hits,$motifs,$seq)=@_;
  my $hash=$motifs->{$seq};
  my @keys=keys %$hash;
  foreach my $key (@keys) {
    $hits->{$key}=1;
  }
}

