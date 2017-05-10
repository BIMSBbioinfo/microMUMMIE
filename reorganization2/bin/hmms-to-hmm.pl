#!/usr/bin/env perl
use strict;

####################################################################
# Converts *.hmms files, which are human-readable and editable, to
# *.hmm files, which are more easily machine-processable.
#
# bmajoros@tigr.org
####################################################################

my $usage="$0 <*.hmms> <*.hmm>";
die "$usage\n" unless @ARGV==2;
my ($hmmsFile,$hmmFile)=@ARGV;

# READ HMMS FILE
my @alphabet=('A','C','G','N','T');
my (%emit,%trans,$numStates);
open(HMMS,$hmmsFile) || die "Can't open file: $hmmsFile\n";
while(<HMMS>)
  {
    if(/(\d+)\s*->\s*(\d+)\s*:\s*(\S+)/)
      {
	my ($fromState,$toState,$P)=($1,$2,$3);
	$trans{$fromState}->{$toState}=$P;
	if($fromState+1>$numStates) {$numStates=$fromState+1}
	if($toState+1>$numStates) {$numStates=$toState+1}
      }
    elsif(/^\s*state\s+(\d+)\s*:\s*(\S.*\S)\s*$/)
      {
	my ($state,$emissions)=($1,$2);
	my @emissions=split/\s+/,$emissions;
	foreach my $emission (@emissions)
	  {
	    $emission=~/(\S)=(\S+)/ || die "Syntax error in line: $_";
	    my ($symbol,$P)=($1,$2);
	    $emit{$state}->{$symbol}=$P;
	  }
      }
  }
close(HMMS);

# OUTPUT ALPHABET
open(HMM,">$hmmFile") || die "Can't create file: $hmmFile\n";
#print HMM "ACGNT\n";

# OUTPUT NUMBER OF STATES
print HMM "$numStates\n$numStates $numStates\n";

# OUTPUT STATE TRANSITION TABLE
for(my $toState=0 ; $toState<$numStates ; ++$toState)
  {
    for(my $fromState=0 ; $fromState<$numStates ; ++$fromState)
      {
	my $P=0+$trans{$fromState}->{$toState};
	print HMM "$P\t";
      }
    print HMM "\n";
  }
print HMM "\n";

# OUTPUT EMISSION TABLE
print HMM "$numStates 5\n";
foreach my $symbol (@alphabet)
  {
    for(my $state=0 ; $state<$numStates ; ++$state)
      {
	my $P=0+$emit{$state}->{$symbol};
	print HMM "$P\t";
      }
    print HMM "\n";
  }

close(HMM);



