#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;
use FastaReader;
use Translation;
use TempFilename;

my $EXTEND_FOR_FOLDING=200; # number of bp to extend upstream
my $MAX_CHROM_NAME_LEN=3;
my $MAX_GENE_LENGTH=1000000;
my $TWO_BIT_FILE="Genome.2bit";
my $DATASETS_SUBDIR="datasets";
my $EXON_FILE="exons.txt";
my $UTRS_FILE="UTRs.txt";
#my $UTRS_FILE="Hi.Expressed.3UTRs.txt";

my $USAGE="get-folding-sequence.pl <outdir>\n";
die $USAGE unless @ARGV==1;
my $BASEDIR=shift @ARGV;
my @libs=@ARGV;

# LOAD THE EXON FILE
my %exons; # indexed by transcript ID
open(IN,$EXON_FILE) || die $EXON_FILE;
while(<IN>) {
  chomp;
  my @fields=split/\s+/,$_;
  next unless @fields==3;
  my ($id,$begin,$end)=@fields;
  --$begin;
  my $rec={begin=>$begin,end=>$end};
  push @{$exons{$id}},$rec;
}
close(IN);

# LOAD THE UTR FILE
my %genes;
my $bad=0;
open(IN,$UTRS_FILE) || die $UTRS_FILE;
while(<IN>) {
  chomp;
  my @fields=split/\t/,$_;
  my ($chr,$begin,$end,$geneID,$strand,$transID)=@fields;
  --$begin;
  next unless length($chr)<=$MAX_CHROM_NAME_LEN;
  my $gene=$genes{$geneID};
  if(!defined($gene)) {
    $gene=$genes{$geneID}=
      {
       chr=>$chr,
       strand=>$strand,
      };
  }
  if($chr ne $gene->{chr} || $strand ne $gene->{strand})
    {
      ++$bad;
      delete $genes{$geneID};
      next;
    }
  my $utr=
    {
     begin=>$begin,
     end=>$end,
     transcriptID=>$transID
    };
  push @{$gene->{UTRs}},$utr;
}
close(IN);

my @geneIDs=keys %genes;
my $numGenes=@geneIDs;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $geneID=$geneIDs[$i];
  my $gene=$genes{$geneID};
  my $chr=$gene->{chr};
  my $strand=$gene->{strand};
  my $groups=$gene->{groups};
  my $UTRs=$gene->{UTRs};
  my $numUTRs=@$UTRs;
  my (%isoforms,%exonIDs);
  my $nextIsoID=1;
  my $geneMin=$UTRs->[0]->{begin};
  my $geneMax=$geneMin;
  for(my $i=0 ; $i<$numUTRs ; ++$i) {
    my $utr=$UTRs->[$i];
    my $utrBegin=$utr->{begin};
    my $utrEnd=$utr->{end};
    if($utrBegin<$geneMin) {$geneMin=$utrBegin}
    if($utrEnd>$geneMax) {$geneMax=$utrEnd}
  }
  if($geneMax-$geneMin>$MAX_GENE_LENGTH) {next}
  if($strand eq "+") { @$UTRs=sort {$a->{begin} <=> $b->{begin}} @$UTRs }
  else { @$UTRs=sort {$b->{begin} <=> $a->{begin}} @$UTRs }
  my %byTranscript;
  for(my $i=0 ; $i<$numUTRs ; ++$i) {
    my $utr=$UTRs->[$i];
    my $transcriptID=$utr->{transcriptID};
    push @{$byTranscript{$transcriptID}},$utr;
  }
  my @transcriptIDs=keys %byTranscript;
  foreach my $transcriptID (@transcriptIDs) {
    my $exons=$exons{$transcriptID};
    if($strand eq "+") {@$exons=sort {$a->{begin} <=> $b->{begin}} @$exons}
    else {@$exons=sort {$b->{begin} <=> $a->{begin}} @$exons}
    my $UTRs=$byTranscript{$transcriptID};
    extend($UTRs,$exons,$strand,$transcriptID) || next;
    my $numUTRs=@$UTRs;
    my $totalSeq;
    for(my $i=0 ; $i<$numUTRs ; ++$i) {
      my $utr=$UTRs->[$i];
      my $utrBegin=$utr->{begin};
      my $utrEnd=$utr->{end};
      my $tempFile=TempFilename::generate();
      my $msg=`twoBitToFa $TWO_BIT_FILE $tempFile -seq=chr$chr -start=$utrBegin -end=$utrEnd -noMask`; # half-open zero-based coordinates are right
      if($msg=~/is not in/) {system("rm -f $tempFile");next}
      if(-z "$tempFile") {system("rm -f $tempFile");next}
      my $reader=new FastaReader("$tempFile");
      my ($extdef,$extseq)=$reader->nextSequence();
      if($strand eq "-")
	{ $extseq=Translation::reverseComplement(\$extseq) }
      $reader->close();
      system("rm -f $tempFile");
      $totalSeq.=$extseq;
    }
    my $outfile="$BASEDIR/$transcriptID.fasta";
    open(OUT,">$outfile") || die $outfile;
    #$seq="\U$seq";
    print OUT ">extendedDNA /geneID=$geneID /transcriptID=$transcriptID\n";
    print OUT "$totalSeq\n";
    close(OUT);
  }
}

#########################################################################
sub extend {
  my ($UTRs,$exons,$strand,$transID)=@_;
  my $firstUTR=$UTRs->[0];
  my $begin=$firstUTR->{begin};
  my $end=$firstUTR->{end};
  my $utrExonIndex;
  my $numExons=@$exons;
  for(my $i=$numExons-1 ; $i>=0 ; --$i) {
    my $exon=$exons->[$i];
    if($begin<$exon->{end} && $exon->{begin}<$end) {
      $utrExonIndex=$i;
      last;
    }
  }
  if(!defined $utrExonIndex) {return 0}
  my $utrExon=$exons->[$utrExonIndex];
  if($strand eq "+") {
    my $newBegin=$begin-$EXTEND_FOR_FOLDING;
    if($newBegin>=$utrExon->{begin}) { $firstUTR->{begin}=$newBegin; return 1 }
    $firstUTR->{begin}=$utrExon->{begin};
    my $deficit=$utrExon->{begin}-$newBegin;
    for(my $i=$utrExonIndex-1 ; $i>=0 ; --$i) {
      my $exon=$exons->[$i];
      my $exonLen=$exon->{end}-$exon->{begin};
      if($exonLen>=$deficit) {
	my $rec={begin=>$exon->{end}-$deficit,
		 end=>$exon->{end},
		 transcriptID=>$transID};
	unshift(@$UTRs,$rec);
	return 1;
      }
      my $rec={begin=>$exon->{begin},
	       end=>$exon->{end},
	       transcriptID=>$transID};
      unshift(@$UTRs,$rec);
      $deficit-=$exonLen;
    }
  }
  else {
    my $newEnd=$end+$EXTEND_FOR_FOLDING;
    if($newEnd<=$utrExon->{end}) { $firstUTR->{end}=$newEnd; return 1 }
    $firstUTR->{end}=$utrExon->{end};

    my $deficit=$newEnd-$utrExon->{end};
    for(my $i=$utrExonIndex-1 ; $i>=0 ; --$i) {
      my $exon=$exons->[$i];
      my $exonLen=$exon->{end}-$exon->{begin};
      if($exonLen>=$deficit) {
	my $rec={begin=>$exon->{begin},
		 end=>$exon->{begin}+$deficit,
		 transcriptID=>$transID};
	unshift(@$UTRs,$rec);
	return 1;
      }
      my $rec={begin=>$exon->{begin},
	       end=>$exon->{end},
	       transcriptID=>$transID};
      unshift(@$UTRs,$rec);
      $deficit-=$exonLen;
    }
  }
  return 0;
}


