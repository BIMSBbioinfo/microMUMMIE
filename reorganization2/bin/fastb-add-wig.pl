#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use ProgramName;

#============================================================
my $name=ProgramName::get();
die "$name <indir-fastb> <indir-wig> <outdir>\n" unless @ARGV==3;
my ($fastbDir,$wigDir,$outdir)=@ARGV;
my @fastbFiles=`ls $fastbDir`;
foreach my $fastb (@fastbFiles) {
  $fastb=~/([^\/]+)\.fastb/ || next;
  my $stem=$1;
  my $fastbFile="$fastbDir/$stem.fastb";
  open(FASTB,$fastbFile) || die "can't open $fastbFile\n";
#print "$fastbFile\n";
  my $outfile="$outdir/$stem.fastb";
  open(OUT,">$outfile") || die "can't write to file $outfile\n";
  my $transcriptID;
  while(<FASTB>) {
    if(/transcriptID=(\S+)/) {$transcriptID=$1}
    print OUT;
  }
  close(FASTB);
  my $wigFile="$wigDir/$transcriptID.fastb";
#die $wigFile unless -e $wigFile;
  my $L=`wc -l $wigFile`;
  my $L1=$L-1;
  #if(!open(WIG,"tail -$L1 $wigFile | add-abscissa.pl | smooth.pl 3 5|")) {
  if(!open(WIG,"tail -$L1 $wigFile | add-abscissa.pl |")) {
    close(OUT);
    system("rm $outfile");
    next;
  }
  my $header=`head -1 $wigFile`;
  chomp $header;
  #while(<WIG>) { print OUT }
  my @wig;
  push @wig,$header;
  while(<WIG>) { if(/\S+\s+(\S+)/){ push @wig,$1 }}
  close(WIG);
  #pop @wig;###
  foreach my $line (@wig) {print OUT "$line\n"}
  close(OUT);
}




