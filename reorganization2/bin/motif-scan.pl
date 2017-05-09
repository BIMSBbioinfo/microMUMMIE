#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;
use FastaReader;
use Fastb;

my $name=ProgramName::get();
die "$name <motifs.fasta> <fastb-dir> <ignore-margin>\n" unless @ARGV==3;
my ($motifFile,$dir,$margin)=@ARGV;

my %motifs;
my $reader=new FastaReader($motifFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/>(\S+)/ || die;
  #$motifs{$sequence}=$1;

  ### miRNA-specific:
  push @{$motifs{substr($sequence,0,8)}},$1;      # 8mer-m1
  push @{$motifs{substr($sequence,0,7)."A"}},$1;  # 8mer-A1
  push @{$motifs{substr($sequence,1,7)}},$1;      # 7mer-m1
  push @{$motifs{substr($sequence,1,6)."A"}},$1;  # 7mer-A1
  push @{$motifs{substr($sequence,0,7)}},$1;      # 7mer-m8
  push @{$motifs{substr($sequence,0,6)}},$1;      # 6mer3-8
  push @{$motifs{substr($sequence,1,6)}},$1;      # 6mer1-7
}

#open(OUTL,">scan-left.txt") || die;
#open(OUTR,">scan-right.txt") || die;
#open(OUTN,">scan-none.txt") || die;
my $fh;

my @files=`ls $dir`;
foreach my $file (@files) {
  $file=~/(.*)\.fastb/ || next;
  my $gene=$1;
  my $fastb=new Fastb("$dir/$file");
  my $dnaTrack=$fastb->getTrackByName("DNA");
  my $dna=$dnaTrack->getData();
  my $numTracks=$fastb->numTracks();

  my $L=length($dna);
  my $first=$margin;
  my $last=$L-8-$margin;
  for(my $i=$first ; $i<=$last ; ++$i) {
    my $substr=substr($dna,$i,8);
    my $hits=$motifs{$substr};
    if($hits && @$hits) {
      #my $c=classifyHit(\$dna,$i,$i+8,$L);
      #if($c<0) {$fh=\*OUTL}
      #elsif($c>0) {$fh=\*OUTR}
      #else {$fh=\*OUTN}
      my $begin=$i+1;
      my $end=$i+8;
      foreach my $hit (@$hits) {
	my $extra=getExtraFields($fastb,$numTracks,$i,$i+8);
	$extra.="seq=$substr;miRNA=$hit;";
	#print $fh "$gene\t$hit\tsite\t$begin\t$end\t8\t+\t.\t$extra\n";
	print "$gene\tnaive\tsite\t$begin\t$end\t8\t+\t.\t$extra\n";
      }
    }
    my $substr=substr($dna,$i,7);
    my $hits=$motifs{$substr};
    if($hits && @$hits) {
      #my $c=classifyHit(\$dna,$i,$i+7,$L);
      #if($c<0) {$fh=\*OUTL}
      #elsif($c>0) {$fh=\*OUTR}
      #else {$fh=\*OUTN}
      my $begin=$i+1;
      my $end=$i+7;
      foreach my $hit (@$hits) {
	my $extra=getExtraFields($fastb,$numTracks,$i,$i+7);
	$extra.="seq=$substr;miRNA=$hit;";
	#print $fh "$gene\t$hit\tsite\t$begin\t$end\t7\t+\t.\t$extra\n";
	print "$gene\tnaive\tsite\t$begin\t$end\t7\t+\t.\t$extra\n";
      }
    }
    my $substr=substr($dna,$i,6);
    my $hits=$motifs{$substr};
    if($hits && @$hits) {
      #my $c=classifyHit(\$dna,$i,$i+6,$L);
      #if($c<0) {$fh=\*OUTL}
      #elsif($c>0) {$fh=\*OUTR}
      #else {$fh=\*OUTN}
      my $begin=$i+1;
      my $end=$i+6;
      foreach my $hit (@$hits) {
	my $extra=getExtraFields($fastb,$numTracks,$i,$i+6);
	$extra.="seq=$substr;miRNA=$hit;";
	#print $fh "$gene\t$hit\tsite\t$begin\t$end\t6\t+\t.\t$extra\n";
	print "$gene\tnaive\tsite\t$begin\t$end\t6\t+\t.\t$extra\n";
      }
    }
  }
}
#close(OUTL);
#close(OUTR);
#close(OUTN);

sub getExtraFields {
  my ($fastb,$numTracks,$begin,$end)=@_;
  my $extra="";
  my $len=$end-$begin;
  for(my $i=0 ; $i<$numTracks ; ++$i) {
    my $track=$fastb->getIthTrack($i);
    next unless $track->isContinuous();
    my $name=$track->getID();
    my $sum=0;
    for(my $j=$begin ; $j<$end ; ++$j) {
      $sum+=$track->getData()->[$j];
    }
    $sum/=$len;
    $extra.="$name=$sum;";
  }
  return $extra;
}

sub isT {
  my ($dna,$pos,$L)=@_;
  return $pos>=0 && $pos<$L && substr($$dna,$pos,1) eq "T";
}

sub classifyHit {
  my ($dna,$begin,$end,$L)=@_;
  for(my $dist=1 ; $dist<10 ; ++$dist) {
    my $leftIsT=isT($dna,$begin-$dist,$L);
    my $rightIsT=isT($dna,$end+$dist-1,$L);
    if($leftIsT && !$rightIsT) {return -1}
    if($rightIsT) {return 1}
  }
  return 0;
}







