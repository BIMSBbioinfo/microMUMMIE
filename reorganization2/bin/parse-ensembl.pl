#!/usr/bin/perl
# Copyright (C)2012 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
use strict;

my $MIN_DIFF=20; # min differences (bp) between distinct transcripts

print STDERR "reading input...\n";
my (%exons,%byGene);
open(IN,"ensembl.txt") || die;
<IN>;
while(<IN>) {
  chomp;
  my @fields=split/,/,$_;
  my ($EnsemblGeneID,$EnsemblTranscriptID,$EnsemblProteinID,$ChromosomeName,
      $GeneStart,$GeneEnd,$TranscriptStart,$TranscriptEnd,$Strand,
      $AssociatedGeneName,$AssociatedGeneDB,$FiveUTRStart,$FiveUTREnd,
      $ThreeUTRStart,$ThreeUTREnd,$CDSLength,$Transcriptcount,$GeneBiotype,
      $EnsemblExonID,$ExonChrStart,$ExonChrEnd,$ConstitutiveExon,
      $ExonRankinTranscript,$phase,$cDNAcodingstart,$cDNAcodingend,
      $Genomiccodingstart,$Genomiccodingend,$CDSStart,$CDSEnd)=@fields;
  if($ExonChrStart>0 && $ExonChrEnd>0) {
    my $rec=
      {
       begin=>$ExonChrStart,
       end=>$ExonChrEnd
      };
    push @{$exons{$EnsemblTranscriptID}},$rec;
  }
  next unless $ThreeUTRStart>0 && $ThreeUTREnd>0;
  my $strand=$Strand>0 ? "+" : "-";
  my $rec=
    {
     chr=>$ChromosomeName,
     utrBegin=>$ThreeUTRStart,
     utrEnd=>$ThreeUTREnd,
     geneName=>$AssociatedGeneName,
     strand=>$strand,
    };
  push @{$byGene{$EnsemblGeneID}->{$EnsemblTranscriptID}},$rec;
}
close(IN);

print STDERR "de-duplicating...\n";
my @geneIDs=keys %byGene;
my $numGenes=@geneIDs;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $geneID=$geneIDs[$i];
  my $gene=$byGene{$geneID}; # hash indexed by transcript ID's
  my @transIDs=keys %$gene;
  my $numTranscripts=@transIDs;
 OUTER:
  for(my $j=0 ; $j<$numTranscripts ; ++$j) {
    my $transcriptID=$transIDs[$j];
    next unless $transcriptID;
    my $transcript=$gene->{$transcriptID}; # array of UTRs
    if(!$transcript) {die $transcriptID}
    my $numUTRs=@$transcript;
    for(my $k=$j+1 ; $k<$numTranscripts ; ++$k) {
      my $otherID=$transIDs[$k];
      next unless $otherID;
      my $otherTrans=$gene->{$otherID};
      my $otherNumUTRs=@$otherTrans;
      if($otherNumUTRs!=$numUTRs) {next}
      if(differences($transcript,$otherTrans)<$MIN_DIFF) {
	$gene->{$transcriptID}=unify($transcript,$otherTrans);
	delete $gene->{$otherID};
	$transIDs[$k]=undef;
	next OUTER;
      }
    }
  }
}
open(OUT,">exons.txt");
my $totalTranscripts=0;
for(my $i=0 ; $i<$numGenes ; ++$i) {
  my $geneID=$geneIDs[$i];
  my $gene=$byGene{$geneID}; # hash indexed by transcript ID's
  my @transIDs=keys %$gene;
  my $numTranscripts=@transIDs;
  $totalTranscripts+=$numTranscripts;
  for(my $j=0 ; $j<$numTranscripts ; ++$j) {
    my $transID=$transIDs[$j];
    my $transcript=$gene->{$transID};
    next unless $transcript;
    my $numExons=@$transcript;
    for(my $k=0 ; $k<$numExons ; ++$k) {
      my $t=$transcript->[$k];
      my $begin=$t->{utrBegin};#-1;
      print "$t->{chr}\t$begin\t$t->{utrEnd}\t$t->{geneName}\t$t->{strand}\t$transID\n";
    }
    my $exons=$exons{$transID};
    foreach my $exon (@$exons) {
      my $begin=$exon->{begin};#-1;
      my $end=$exon->{end};
      print OUT "$transID\t$begin\t$end\n";
    }
  }
}
close(OUT);
#print "total: $totalTranscripts\n";

################################################################
sub transcriptExtent
  {
    my ($transcript)=@_;
    my $n=@$transcript;
    my $exon=$transcript->[0];
    my $begin=$exon->{utrBegin};
    my $end=$exon->{utrEnd};
    for(my $i=1 ; $i<$n ; ++$i) {
      my $exon=$transcript->[$i];
      if($exon->{utrBegin}<$begin) {$begin=$exon->{utrBegin}}
      if($exon->{utrEnd}>$end) {$end=$exon->{utrEnd}}
    }
    return ($begin,$end);
  }
sub differences
  {
    my ($this,$that)=@_;
    my $diff=0;
    my $numExons=@$this;
    for(my $i=0 ; $i<$numExons ; ++$i) {
      my $thisExon=$this->[$i];
      my $thatExon=$that->[$i];
      $diff+=abs($thisExon->{utrBegin}-$thatExon->{utrBegin});
      $diff+=abs($thisExon->{utrEnd}-$thatExon->{utrEnd});
    }
    return $diff;
  }
sub unify
  {
    my ($this,$that)=@_;
    my $diff=0;
    my $numExons=@$this;
    for(my $i=0 ; $i<$numExons ; ++$i) {
      my $thisExon=$this->[$i];
      my $thatExon=$that->[$i];
      my $begin=min($thisExon->{utrBegin},$thatExon->{utrBegin});
      my $end=max($thisExon->{utrEnd},$thatExon->{utrEnd});
      $thisExon->{utrBegin}=$thatExon->{utrBegin}=$begin;
      $thisExon->{utrEnd}=$thatExon->{utrEnd}=$end;
    }
    return $diff;
  }
sub min {
  my ($a,$b)=@_;
  return $a<$b ? $a : $b;
}
sub max {
  my ($a,$b)=@_;
  return $a>$b ? $a : $b;
}
sub overlap
  {
    my ($this,$next)=@_;
    return $this->[0]<$next->[1] && $next->[0]<$this->[1];
  }




