#!/usr/bin/env perl
#073112 -- add replace group num (in cs); ****and join after changing (wasn't like that before!!!)
#070512 -- add replace pcentile (b/c we run separately!)
#052412  -- combined cs and pct
use strict;
use warnings;

my $csF = $ARGV[0];
my $pctF = $ARGV[1];
my $outF = $ARGV[2]; #combined res
my $selsp = $ARGV[3];  #9606 hg
#my $billfolder = $ARGV[4];

my ($line,$pctH,$tcH);

open IN,$pctF or die;
$line=<IN>;

while ($line=<IN>)
{
  chomp $line;
 
  my @a=split(/\t/,$line);

  if ($a[2] ne $selsp) { next; }
  my $id = $a[0]."\t".$a[1]."\t".$a[5]."\t".$a[6]."\t".$a[8]; 
  $pctH->{$id}{"bl"} = $a[11];
  $pctH->{$id}{"pct"} = $a[12];
  if(defined $a[13]) { $pctH->{$id}{"iscs"} = "csv"; } else { $pctH->{$id}{"iscs"} = "ncsv"; } 

  $tcH->{$a[0]} = 1;

#  print $pctH->{$id}{"bl"},"=",$pctH->{$id}{"pct"},"==",$pctH->{$id}{"iscs"},"\n"; 
}

close IN;

#==> tgscn_blpct.out <==
#Gene_ID miRNA_family_ID species_ID      MSA_start       MSA_end UTR_start       UTR_end Group_num       Site_type       miRNA in this species   Group_type      Branch length   Pct     Conserved
#ENST00000327835 ebv-mir-bart1-3p        30608   20366   20372   12217   12223   1       7mer-1a         7mer-1a 0.08523 NA
###################################
open IN,$csF or die;
$line=<IN>;	chomp $line;

open OUT,">".$outF or die;
print OUT $line,"\tBranch length\tPct\tConserved\n";
#foreach my $k (keys %$tcH) 
#{ open B,">".$billfolder."/".$k.".gff" or die; print B $line,"\tBranch length\tPct\tConserved\n"; }

while ($line=<IN>)
{
  chomp $line;
  my @a=split(/\t/,$line);

  if ($a[1] ne $selsp) { next; }
  my $id = $a[0]."\t".$a[16]."\t".$a[4]."\t".$a[5]."\t".$a[3];
  $a[12] = "NA";  #070512 replace pcentle w/ NA b/c sep run
  $a[17] = "NA";  #073112 replace group number to NA as well!!!!!1

  #my $mibid = $a[2];  #diff from a[16] which is mir fam (example just used one mi per fam)

  if (defined $pctH->{$id})
  { 
    print OUT join("\t",@a),"\t",$pctH->{$id}{"bl"},"\t",$pctH->{$id}{"pct"},"\t",$pctH->{$id}{"iscs"},"\n"; 
    ###073112; wasn't joining line before today (so pct was never replaced**)

 #  open B,">>".$billfolder."/".$a[0].".gff" or die; 
  # print B $line,"\t",$pctH->{$id}{"bl"},"\t",$pctH->{$id}{"pct"},"\t",$pctH->{$id}{"iscs"},"\n"; 
   #close B;
  }

  else { print "no bl-pct for $id \n"; }
}

close IN;
close OUT;

#==> tgscn_cs.out <==
#Gene ID Species ID      Mirbase ID      Site Type       UTR start       UTR end 3' pairing contribution local AU contribution   position contribution   TA contribution SPS contribution     context+ score   context+ score percentile       UTR region      UTR-miRNA pairing       mature miRNA sequence   miRNA family    Group #
#ENST00000327835 9606    ebv-mir-bart1-3p        7mer-1a 12366   12372   0.004   -0.028  0.009   -0.036  -0.066  -0.191  65        AUUUUUAAAUACACUGGUGCUAA                        ||||||        CUGUAUCACCUAUCGCCACGAU         ebv-mir-bart1-3p        1


