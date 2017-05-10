#!/usr/bin/env perl
use strict;
use FastaReader;
use Translation;
use TempFilename;

my $PAD=15;

die "par2fastb.pl <distribution.csv> <genome.2bit> <outdir>\n" unless @ARGV==3;
my ($distrFile,$genomeFile,$outdir)=@ARGV;

open(IN,$distrFile) || die "can't open file: $distrFile\n";
while(<IN>) {
  chomp;
  my @fields=split/,/,$_;
  my ($chr,$strand,$begin,$end,$groupID,$signalType)=@fields;
  next unless $signalType eq "Signal";
  for(my $i=0 ; $i<6 ; ++$i) {shift @fields}
  my $data=[];
  @$data=@fields;
  my $n=@$data;
  my $uniform=1;
  my $max=0;
  for(my $i=0 ; $i<$n ; ++$i)
    { if($data->[$i]=~/nan/i) { $data->[$i]=0 } }
  for(my $i=0 ; $i<$n ; ++$i) { if($data->[$i]>$max) {$max=$data->[$i]} }
  for(my $i=0 ; $i<$n ; ++$i) { if($data->[$i]!=$max) {$uniform=0} }
  next if $uniform;
  if($max<=0) {die @$data}
  for(my $i=0 ; $i<$n ; ++$i) { $data->[$i]/=$max }
  $begin-=$PAD;
  $end+=$PAD+1;
  my $outfile="$outdir/$chr.$begin$strand$end.fastb";
  my $tempFile=TempFilename::generate();
  #die $tempFile;
  system("rm $tempFile") if -e $tempFile;
  my $msg=`twoBitToFa $genomeFile $tempFile -seq=$chr -start=$begin -end=$end -noMask`;
  if($msg=~/is not in/) {next}
  if(-z $tempFile) {system("rm $tempFile"); next}
  my $reader=new FastaReader($tempFile);
  my ($def,$seq)=$reader->nextSequence();
  if($strand eq "-") { $seq=Translation::reverseComplement(\$seq) }
  $reader->close();
  system("rm $tempFile");
  open(OUT,">$outfile") || die $outfile;
  print OUT ">DNA /chr=$chr /begin=$begin /end=$end /strand=$strand\n$seq\n";
  my @seq;
  my $dataLen=@$data;
  my $newLen=$dataLen+2*$PAD;
  for(my $i=0 ; $i<$newLen ; ++$i) { $seq[$i]=0 }
  for(my $i=0 ; $i<$dataLen ; ++$i ) { $seq[$i+$PAD]=$data->[$i] }
  $dataLen=$newLen;
  print OUT "\%AGO-Signal /length=$dataLen\n";
  if($strand eq "-") {
    my @tmp;
    for(my $i=0 ; $i<$dataLen ; ++$i) {unshift @tmp,$seq[$i]}
    @seq=@tmp;
  }
  for(my $i=0 ; $i<$dataLen ; ++$i) {print OUT "$seq[$i]\n"}
}
close(OUT);
close(IN);
