#!/usr/bin/env perl
use strict;
use ProgramName;
use TempFilename;
use FastaReader;

my $name=ProgramName::get();
die "$name <chunks> <seeds.fasta> <bulge-type:1/2/3> <signal-variance>\n"
  unless @ARGV==4;
my ($chunksDir,$seedsFile,$bulgeType,$var)=@ARGV;

system("date");

my $tmpdir=TempFilename::generate();
system("mkdir $tmpdir");

my $hash=FastaReader::readAllAndKeepDefs($seedsFile); # hash : id->[def,seq]
open(OUT,">$tmpdir/seeds.fasta") || 
  die "can't create file: $tmpdir/seeds.fasta";
my @keys=keys %$hash;
foreach my $key (@keys) {
  my $pair=$hash->{$key};
  my ($def,$seq)=@$pair;
  my $six=substr($seq,1,6);
  my $firstTwo=substr($six,0,2);
  my $firstThree=substr($six,0,3);
  my $firstFour=substr($six,0,4);
  my $lastFour=substr($six,2,4);
  my $lastThree=substr($six,3,3);
  my $lastTwo=substr($six,4,2);
  my $second=substr($six,1,1);
  my $third=substr($six,2,1);
  my $fourth=substr($six,3,1);
  my $form1="$firstTwo$second$lastFour";
  my $form2="$firstThree$third$lastThree";
  my $form3="$firstFour$fourth$lastTwo";
  my $bulgeSeed;
  if($bulgeType==1)    { $bulgeSeed=$form1 }
  elsif($bulgeType==2) { $bulgeSeed=$form2 }
  elsif($bulgeType==3) { $bulgeSeed=$form3 }
  else { die "bulge type must be 1, 2, or 3" }
  print OUT "$def$bulgeSeed\n";
}
close(OUT);

system("bulge/newrunbulgeprediction.pl $chunksDir $tmpdir/seeds.fasta $tmpdir $bulgeType $var");

system("rm -r $tmpdir");

my $bulgeTypeRN;
for(my $i=0 ; $i<$bulgeType ; ++$i) { $bulgeTypeRN.="I" }

my (%outputLines,%outputNames);
my $filename="bulge-predictions-var$var.gff";
open(IN,$filename) || die $filename;
while(<IN>) {
  chomp;
  my @fields=split/\t/,$_;
  my $extra=$fields[8];
  my (%names,$newExtra);
  while($extra=~/miRNA=([^;]+);/) { 
    $names{$1}=1;
    $extra=~s/miRNA=([^;]+);//;
  }
  my @names=keys %names;
  #foreach my $name (@names) { $newExtra.="miRNA=$name;" }
  if($extra=~/seq=([^;]+);/) { $newExtra.="seq=$1;" }
  my $firstName=$names[0];
  $fields[1]=$firstName;
  $fields[2]="bulge$bulgeTypeRN";
  $fields[8]=$newExtra;
  my $line=join("\t",@fields);
  my $substrate=$fields[0];
  my $begin=$fields[3];
  my $end=$fields[4];
  my $key="$substrate $begin $end";
  if(!defined($outputLines{$key})) { $outputLines{$key}=$line }
  foreach my $name (@names) { $outputNames{$key}->{$name}=1 }
}
close(IN);

open(OUT,">$filename") || die $filename;
my @keys=keys %outputLines;
foreach my $key (@keys) {
  my $line=$outputLines{$key};
  print OUT "$line";
  my $namesHash=$outputNames{$key};
  my @names=keys %$namesHash;
  foreach my $name (@names) { print OUT "miRNA=$name;" }
  print OUT "\n";
}
close(OUT);

system("date");

