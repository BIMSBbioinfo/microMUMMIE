#!/usr/bin/env perl
#072012 -- write from scratch
	#select 23 species, map species to spid, concat chunks; 
	#add .. as missing; write in tgscn format

use strict;	use warnings;

my $orimafF = $ARGV[0];
my $wantedspidF = $ARGV[1]; #mappedspid23w.txt  for tgscn*
my $wantedenstF = $ARGV[2]; #../wtconslongest/wtclgst_enst.txt
my $outF = $ARGV[3];

my ($line,$enstH,$spidH,$H);

#-------get wanted spid map----------------------------------------------------- 
open IN,$wantedspidF or die;

while ($line = <IN>)
{
  chomp $line;	if ($line eq "") { next; }
  my @a = split ("\t",$line);
 
  $spidH->{lc($a[0])} = $a[1];	
}

close IN;
#hg19    9606

#-------get wanted enst--------------------------------------------------------- 
open IN,$wantedenstF or die;

while ($line = <IN>)
{
  chomp $line;	if ($line eq "") { next; }
  $enstH->{$line} = 1;	
}

close IN;

#ENST00000443805

#-------get maf and parse ------------------------------------------------------ 
open IN,$orimafF or die;
my ($spid,$enst,$pf,$sf);

while ($line = <IN>)
{
  chomp $line;	if ($line eq "") { next; }
  
  if ($line =~ m/^#/) { next; }
 
  if ($line =~ m/^a/)
  {
    $line = <IN>;  chomp $line;
    my @a = split(" ",$line); 
    my @a2 =split("_",$a[1]);  $sf = pop @a2; $pf=join("_",@a2);
    $enst = $a[1];
    
    print $enst,"\n";
    
    if (!defined $enstH->{$pf}) {next; }
    
    $spid = "9606";
    $H->{$pf}{$sf}{$spid} = $a[6];
    
    next;
  }
  
  if ($line =~ m/^s/)  #not main anchor (enst)
  {
    if (!defined $enstH->{$pf}) { next; }

    my @a = split(" ",$line);
    if (!defined $spidH->{lc($a[1])} ) { next; }
    
    $spid = $spidH->{lc($a[1])};
    $H->{$pf}{$sf}{$spid} = $a[6];
  }
}

close IN;

##maf version=1 scoring=zero
#a score=0.000000
#s ENST00000000233_1 0 406 + 406 C--CAGC-C----------A-GGGGCAGGCCCCTGATGCCCGGAAGCTCCTGCG-TG---CA------T---CCCCGGGATGACCAGACTCCC--GGACTCCTCAGGC---------AGTGCCCTTTCCTC--C----CA-CTTTTCCTCCCC------------------------------------CAT--AG-----CCACAGGC----C----TCTGCTC--------C------------TG---CTCCTGCCT-GCATGTTCTCTC-TGT-TGTTGGAGCC--TG-GAG-----------


#------
my $cbnseqH;
foreach my $pfk (sort keys %$H)
{
  foreach my $spk ( sort keys %$spidH )
  {
    my $spid = $spidH->{$spk};
    $cbnseqH->{$pfk}{$spid} = "";
  }
}


foreach my $pfk (sort keys %$H)
{
  foreach my $sfk (sort {$a<=>$b} keys %{$H->{$pfk}})
  {
    my $anclen = length ($H->{$pfk}{$sfk}{"9606"});
    $cbnseqH->{$pfk}{"9606"} .= $H->{$pfk}{$sfk}{"9606"};
     
    foreach my $spk ( sort keys %$spidH )
    {
      my $spid = $spidH->{$spk};
      
      if ($spid eq "9606") { next; }
    
      if (!defined $H->{$pfk}{$sfk}{$spid}) 
      {
        for (my $i=0; $i<$anclen; $i++) { $cbnseqH->{$pfk}{$spid} .= "."; }
      }
      
      else
      {
        if (length ($H->{$pfk}{$sfk}{$spid}) ne $anclen) 
        { print "diff len in alm: $pfk $sfk $spid\n"; die; }
        
        $cbnseqH->{$pfk}{$spid} .= $H->{$pfk}{$sfk}{$spid};
      }
    }
  }
}

open OUT,">".$outF or die;

foreach my $k1 (sort keys %$cbnseqH)
{
  print OUT $k1,"\t9606\t",$cbnseqH->{$k1}{"9606"},"\n";

  foreach my $k2 (sort keys %{$cbnseqH->{$k1}})
  {
    if ($k2 eq "9606") { next; }
    print OUT $k1,"\t",$k2,"\t",$cbnseqH->{$k1}{$k2},"\n";
  }
}

close OUT;

#name\t9606\tseq
