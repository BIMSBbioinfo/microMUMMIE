#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;
use FastaReader;

my $INCLUDE_8MERS=1;
my $INCLUDE_7MERS=1;
my $INCLUDE_6MERS=1;

my $name=ProgramName::get();
die "$name <in.gff> <seeds.fasta>\n" unless @ARGV==2;
my ($gffFile,$seedsFile)=@ARGV;

my (%motifs1_8,%motifs1_7,%motifs2_8,%motifs2_7,%motifs1_6);
my $reader=new FastaReader($seedsFile);
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  $defline=~/>(\S+)/ || die;
  $motifs1_8{substr($sequence,0,8)}->{$1}=1;      # 8mer-m1
  $motifs1_8{substr($sequence,0,7)."A"}->{$1}=1;  # 8mer-A1
  $motifs2_8{substr($sequence,1,7)}->{$1}=1;      # 7mer-m1
  $motifs2_8{substr($sequence,1,6)."A"}->{$1}=1;  # 7mer-A1
  $motifs1_7{substr($sequence,0,7)}->{$1}=1;      # 7mer-m8
  $motifs1_6{substr($sequence,0,6)}->{$1}=1;      # 6mer3-8
  $motifs2_7{substr($sequence,1,6)}->{$1}=1;      # 6mer2-7
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
    if($field=~/path=([^;]+)/) { @path=split/,/,$1 }
  }
  #die unless @path>0;
  my %hits;
  if(@path==0) {
    if($INCLUDE_8MERS)  { addHits(\%hits,\%motifs1_8,substr($seq,0,8)) }
    if($INCLUDE_7MERS)  { addHits(\%hits,\%motifs1_7,substr($seq,0,7)) }
    if($INCLUDE_7MERS)  { addHits(\%hits,\%motifs2_8,substr($seq,1,7)) }
    if($INCLUDE_6MERS)  { addHits(\%hits,\%motifs2_7,substr($seq,1,6)) }
    if($INCLUDE_6MERS)  { addHits(\%hits,\%motifs1_6,substr($seq,0,6)) }
  }
  else {
    my $startState=$path[0];
    my ($newSeq,$type);
    if($startState==5) {
      $newSeq=substr($seq,0,8);
      addHits(\%hits,\%motifs1_8,$newSeq);
      if($path[7]==12) { $type="8mer-m1" } else { $type="8mer-A1" }
    }
    elsif($startState==14) {
      $newSeq=substr($seq,1,7);
      addHits(\%hits,\%motifs2_8,$newSeq);
      $site->setBegin($site->getBegin()+1);
      if($path[7]==21) { $type="7mer-m1" } else { $type="7mer-A1" }
    }
    elsif($startState==22) {
      $newSeq=substr($seq,0,7);
      addHits(\%hits,\%motifs1_7,$newSeq);
      $site->setEnd($site->getEnd()-1);
      $type="7mer-m8";
    }
    elsif($startState==30) {
      $newSeq=substr($seq,1,6);
      addHits(\%hits,\%motifs2_7,$newSeq);
      $site->setBegin($site->getBegin()+1);
      $site->setEnd($site->getEnd()-1);
      $type="6mer2-7";
    }
    elsif($startState==38) {
      $newSeq=substr($seq,0,6);
      addHits(\%hits,\%motifs1_6,$newSeq);
      $site->setEnd($site->getEnd()-2);
      $type="6mer3-8";
    }
    else { my $path=join(',',@path); die $path }
    $extra->[$seqField]="seq=$newSeq;";
    push @$extra,"type=$type;";
  }
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

