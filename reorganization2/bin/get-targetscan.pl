#!/usr/bin/perl
use strict;
use ProgramName;

#my $N_MEANS_NONCONSERVED=0;

my $name=ProgramName::get();
die "$name <rawfile> <screenfile> <cons:Y/N> <transcript-ids.txt> <nonconserved:0/1>\n" unless @ARGV==5;
my ($RAWFILE,$MIRNA_SCREEN_FILE,$consYN,$TRANSCRIPT_ID_FILE,
   $N_MEANS_NONCONSERVED)=@ARGV;

my $REQUIRE_CONSERVATION=($consYN eq "Y" ? 1 : 0);
#my $MIRNA_SCREEN_FILE="15percent/mirbase/top50.fasta";
#my $RAWFILE="15percent/targetscan/shuffled.txt.gz";
#my $RAWFILE="targetscan/sanity/shuffled.txt.gz";
#my $RAWFILE="targetscan/sanity/unshuffled.txt";
#my $RAWFILE="targetscan/ShuffledAug6.txt.gz";
#my $RAWFILE="targetscan/unshuffled.txt";

my %keepMiRNAs;
open(IN,$MIRNA_SCREEN_FILE) || die $MIRNA_SCREEN_FILE;
while(<IN>) {
  chomp;
  if(/>(\S+)/) {$keepMiRNAs{$1}=1}
}
close(IN);

my %ids;
open(IN,$TRANSCRIPT_ID_FILE) || die $TRANSCRIPT_ID_FILE;
while(<IN>) {
  chomp;
  if(/(\S+)\s+(\S+)/) { $ids{$2}=$1 }
}
close(IN);

#my $RAWFILE="/nfs/igsp/labs/ohlerlab/nicky/dparclip/mvHMM/othertool/tgscn/WT-23walm-1mifam/combined9606_cspct.out";
#my $RAWFILE="/nfs/igsp/labs/ohlerlab/nicky/dparclip/mvHMM/othertool/tgscn/shWT-23walm-1mifam/combined96ombined9606_cspct.out";
#my $RAWFILE="targetscan/utr.txt";
if($RAWFILE=~/\.gz$/) { open(IN,"cat $RAWFILE|gunzip|") || die }
else { open(IN,$RAWFILE) || die }
<IN>;
while(<IN>) {
  chomp;
  my @fields=split/\t/,$_;
  my ($GeneID,$SpeciesID,$MirbaseID,$SiteType,$UTRstart,$UTRend,
      $threePrimePairingContribution,$localAUcontribution,
      $positionContribution,$TAcontribution,$SPScontribution,
      $contextScore,
      $contextScorePercentile,$UTRregion,$UTRmiRNApairing,$matureMiRNAsequence,$miRNAfamily,$GroupNum,
      $BranchLength,
      $Pct,$Conserved)=@fields;
  #print "branch=$BranchLength  pct=$Pct   C=$Conserved\n";
  my $numFields=@fields;
  #next if $Pct eq "NA";
  if($REQUIRE_CONSERVATION) {next unless $Conserved eq "csv"}
  if($N_MEANS_NONCONSERVED && !$REQUIRE_CONSERVATION)
    {next unless $Conserved eq "ncsv"}
  my $lookupID=$MirbaseID;
  $lookupID=~s/_sh_\d+//g;
  next unless $keepMiRNAs{$lookupID};
  my $id=$fields[0];
  my $begin=$fields[4];
  my $end=$fields[5];
  my $score=$fields[12];
  if($id=~/:(\S+)/) {$id=$1}
  if(!defined($ids{$id})) {next} 
  $id=$ids{$id};
  $score=$contextScore;
  print "$id\tTargetScanS\tsite\t$begin\t$end\t$score\t+\t.\tmiRNA=$MirbaseID;$Conserved;score3p=$threePrimePairingContribution;scoreAU=$localAUcontribution;scorePos=$positionContribution;scoreAbund=$TAcontribution;scoreSPS=$SPScontribution;branchLen=$BranchLength;\n";
}
close(IN);

