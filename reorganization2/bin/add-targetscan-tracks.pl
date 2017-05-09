#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <targetscan.gff> <indir> <outdir> <score-type>\n" unless @ARGV==4;
my ($targetscanGFF,$indir,$outdir,$scoreType)=@ARGV;

my ($currentSubstrate,@records);
open(IN,$targetscanGFF) || die;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my ($substrate,$ts,$site,$begin,$end,$score,$strand,$frame,$extra)=@fields;
  next unless $score=~/\d/;
  if(!defined($currentSubstrate)) {$currentSubstrate=$substrate}
  if($substrate ne $currentSubstrate) {
    process($currentSubstrate,\@records);
    @records=();
    $currentSubstrate=$substrate;
  }
  if($scoreType eq "score3p"    ||
     $scoreType eq "scoreAU"    ||
     $scoreType eq "scorePos"   ||
     $scoreType eq "scoreAbund" ||
     $scoreType eq "scoreSPS"   ||
     $scoreType eq "branchLen")
    {$score=extractExtra($extra,$scoreType)}
  push @records,{begin=>$begin,end=>$end,score=>$score};
}
close(IN);
process($currentSubstrate,\@records);

@records=();
my @infiles=`ls $indir`;
foreach my $file (@infiles) {
  chomp $file;
  next if -e "$outdir/$file";
  $file=~/(\S+)\.fastb/ || die;
  my $substrate=$1;
  process($substrate,\@records);
}

sub extractExtra {
  my ($extra,$field)=@_;
  $extra=~/$field=([^;]+);/ || die $extra;
  return $1;
}
sub process {
  my ($substrate,$records)=@_;
  my $infile="$indir/$substrate.fastb";
  return $infile unless -e $infile; 
  my $fastb=new Fastb($infile);
  my $L=$fastb->getLength();
  my $outfile="$outdir/$substrate.fastb";
  my $n=@$records;
  my @sums;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $record=$records->[$i];
    my $begin=$record->{begin}-1;
    my $end=$record->{end};
    my $score=$record->{score};
    for(my $pos=$begin ; $pos<$end ; ++$pos) {
      push @{$sums[$pos]},$score;
    }
  }
  open(OUT,">$outfile") || die;
  open(IN2,$infile) || die;
  while(<IN2>) {print OUT}
  close(IN2);
  print OUT "\%targetscan /length=$L\n";
  for(my $pos=0 ; $pos<$L ; ++$pos) {
    my $array=$sums[$pos];
    my $mean=$array ? mean($array) : 0;
    print OUT "$mean\n";
  }
  close(OUT);
}

sub mean {
  my ($array)=@_;
  my $sum=0;
  my $n=@$array;
  for(my $i=0 ; $i<$n ; ++$i) { $sum+=$array->[$i] }
  $sum/=$n;
  return $sum;
}



