#!/usr/bin/env perl
#Program written to print genomic coordinates from UTR file in the output gff file September 27 2013 Samta Malhotra
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff> <UTRs.txt> <directory fastb files> <out.gff>\n" unless @ARGV==4;
my ($infile,$UTRS_FILE,$dir,$outfile)=@ARGV;
my @data;
my @ensemblid;
open(UTR,$UTRS_FILE) || die "can't open file: $UTRS_FILE\n";
while (<UTR>) {
push @data, $_;
}
my $len=scalar(@data);
my @files=`ls $dir`;
my $n=@files;
open(IN,$infile) || die "can't open file: $infile\n";
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
while(<IN>) {
chomp;
my @fields=split(/\s+/,$_);
my $substrate=$fields[0];
my @value = split(/_/,$substrate);
my $GeneID = $value[0];
my $begin=$fields[3];
my $end=$fields[4];
my $score=$fields[5];
my $strand=$fields[6];
my $dot=$fields[7];
my $extra=$fields[8];
for(my $j=0 ; $j<$n ; ++$j) {
  my $file=$files[$j];
  if ($file =~ m/(^$substrate).fastb/) {
  my $name=$substrate;
  my $sscore = $score;
  open(FILE,"$dir/$substrate.fastb") || die "can't open $dir/$substrate.fastb\n";
  while(<FILE>) {
    
    if($_ =~ m/transcriptID=(ENST\d+)/) 
    {
      my $transcriptid = $1;
      for(my $i=0 ; $i<$len ; ++$i){
      my $line=$data[$i];
      #my $ensemblid = ($line =~ /(ENST\d+)/g);
     # if ($line =~ m/(ENST\d+) /) { my $id = $1;}
     # print my $id;
      my @row=split(/\s+/,$line);
      my $chr = $row[0];
      my $utrstart = $row[1];
      my $utrend = $row[2];
      my $sstrand = $row[4];
      my ($chunkBegin,$chunkEnd)=($1,$2);
      if($sstrand eq "+") {
         $chunkBegin=$utrstart+$begin;
         $chunkEnd=$utrstart+$end;
       }
       else {
       $mystartvalue = $utrend-$begin;
       $myendvalue=$utrend-$end;
       $chunkBegin=$myendvalue;
       $chunkEnd=$mystartvalue;
       }
      my $id = $row[5];
      if ($transcriptid eq $id) {
      my @gene = split(/_/,$name);
      print OUT "chr$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$GeneID;transcriptid=$id;utr_start=$utrstart;utr_end=$utrend;prediction_start=$begin;prediction_end=$end;$extra\n"; last
     }
     else {$line=$data[$i+1];}
    }
   }
  }
  close(FILE);
}
else {$file=$files[$j+1];}
}
#print "$gene\t$id\t$transtart\t$tranend\n";
}
close (IN);
close (OUT);
close (UTR);
