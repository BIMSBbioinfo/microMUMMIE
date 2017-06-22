#!/usr/bin/env perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Changed by Samta 03/03/2014 in order to include the start and end in the description for using it to convert into Genomic coordinates
use strict;
use List::Util qw( min max );

die "assemble-transcripts.pl <indir> <outdir>\n" unless @ARGV==2;
my ($indir,$outdir)=@ARGV;

my %transcripts;
my @files=`ls $indir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  $file=~/(\S+)_(\d+)_(\d+)\.fastb/ || die "can't parse filename: $file\n";
  my ($gene,$transcript,$exon)=($1,$2,$3);
  my $base="$gene\_$transcript";
  $transcripts{$base}++;
}
my @transcripts=keys %transcripts;
my $n=@transcripts;
for(my $i=0 ; $i<$n ;++$i) {
  my $base=$transcripts[$i];
  my $numExons=$transcripts{$base};
  my $outfile="$outdir/$base.fastb\n";
  my %tracks;
  my ($currentTrack,%trackTypes,%deflines,%basebegin,%baseend,@beginall,@endall);
  for(my $j=0 ; $j<$numExons ; ++$j) {
    my $infile="$indir/$base\_$j.fastb";
    my $idtoput = join("_",$base,$j);
    open(IN,$infile) || die "can't open file $infile\n";
    while(<IN>) {
      if(/^([>%])(\S+)(.*)/) {
	my ($type,$id,$rest)=($1,$2,$3);
	$currentTrack=$id;
	$trackTypes{$id}=$type;
	$deflines{$id}=$_;
       #if ($rest =~ m/transcriptID=ENST(\d+)\s+\/chr=(\d+)\s+\/begin=(\d+)\s+\/end=(\d+)/){
          if($rest =~ m/transcriptID=(\D+\d+)/)
      {
           my $trans = $1;
          if ($rest =~ m/chr=(.)/){
          my $chrom = $1;
         if ($rest =~ m/begin=(\d+)/){
          my $begin = $1;
        if ($rest =~ m/end=(\d+)/){
          my $end = $1;
          $basebegin{$idtoput}=$begin;
          $baseend{$idtoput}=$end;
      }
      }
      }
      }


      }
      else {
	$tracks{$currentTrack}.=$_;
      }

    }
    close(IN);

  }

  while ( my ($key, $value) = each(%basebegin) ) {
        #print "$key => $value\n";


    }
    my ($value,$valuend);
    for(my $j=0 ; $j<$numExons ; ++$j) {
    my $idtoput = join("_",$base,$j);
    my @begins = keys %basebegin;
    my @ends = keys %baseend;
    foreach my $keybe (@begins){
    if ($keybe = $idtoput) {
   $value = $basebegin{$keybe};
  }
  }
  push (@beginall,$value);
  foreach my $keyen (@ends){
    if ($keyen = $idtoput) {
   $valuend = $baseend{$keyen};
  }
  }
  push (@endall,$valuend);
  }

my $minstart = min @beginall;
my $minend = max @endall;

  open(OUT,">$outfile") || die "can't write to file $outfile\n";
  my @tracks=keys %tracks;
  foreach my $trackID (@tracks) {
    my $data=$tracks{$trackID};
    my $defline=$deflines{$trackID};
    if($defline=~/^%/) {
      $defline=~/^(.*\s+\/length=)\d+\s*$/ || die "can't parse: $defline";
      my $firstPart=$1;
      my $L=($data=~s/\n/\n/g);
      $defline="$firstPart$L\n";
    }
    else {
      $defline=~s/\s+\/begin=\d+\s+\/end=\d+//g || die "can't parse: $defline";
      chomp $defline;
      for (my $s=0;$s<$numExons;$s++)
      {
      my $start = $beginall[$s];
      my $end = $endall[$s];
      $defline.=" /begin_$s=$start\t/end_$s=$end  ";
      }
      $defline.=" /numexons=$numExons\n";
    }
    print OUT "$defline$data";
  }
  close(OUT);
}
