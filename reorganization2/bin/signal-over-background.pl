#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <indir> <outdir>\n" unless @ARGV==2;
my ($indir,$outdir)=@ARGV;

my @files=`ls $indir`;
foreach my $file (@files) {
  $file=~/([^\/]+)\.fastb/ || die $file;
  my $stem=$1;
  my $infile="$indir/$stem.fastb";
  open(IN,$infile) || die "can't open file $infile\n";
  my (%tracks,$seq,$signal,$background);
  while(<IN>) {
    if(/^[%>](\S+).*/) {
      my $id=$1;
      my $rec=$tracks{$id}=
	{
	 def=>$_,
	 id=>$id,
	 seq=>[]
	};
      $tracks{$id}=$rec;
      $seq=$rec->{seq};
      if($id=~/[Ss]ignal/) {$signal=$seq}
      elsif($id=~/[Bb]ackground/) {$background=$seq}
    }
    elsif(/(\S+)/) {push @$seq,$1}
  }
  close(IN);
  my $outfile="$outdir/$stem.fastb";
  open(OUT,">$outfile") || die "can't write to file $outfile\n";
  my $L=@$signal;
  for(my $i=0 ; $i<$L ; ++$i) { 
    my $bg=$background->[$i];
    if($bg!=0) { $signal->[$i]/=$bg }
  }
  my @tracks=keys %tracks;
  my $N=@tracks;
  for(my $i=0 ; $i<$N ; ++$i) {
    my $trackID=$tracks[$i];
    next if($trackID=~/[Bb]ackground/);
    my $rec=$tracks{$trackID};
    my $def=$rec->{def};
    my $seq=$rec->{seq};
    print OUT "$def";
    my $n=@$seq;
    for(my $i=0 ; $i<$n ; ++$i) {
      my $line=$seq->[$i];
      print OUT "$line\n";
    }
  }
  close(OUT);
}
