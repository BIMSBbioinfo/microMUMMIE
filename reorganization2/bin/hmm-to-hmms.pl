#!/usr/bin/env perl
use strict;

####################################################################
# Converts *.hmm files to *.hmms format, which is more user-readable
# and editable.
#
# bmajoros@tigr.org
####################################################################

my $usage="$0 <*.hmm> <*.hmms>";
die "$usage\n" unless @ARGV==2;
my ($hmmFile,$hmmsFile)=@ARGV;

open(HMM,$hmmFile) || die "Can't open file: $hmmFile\n";
open(HMMS,">$hmmsFile") || die "Can't create file: $hmmsFile\n";

my $alphabet="ACGNT";
my @alphabet=split//,$alphabet;

# READ NUMBER OF STATES
my $numStates=0+<HMM>;

# READ STATE TRANSITION TABLE
my %trans;
<HMM>;
for(my $toState=0 ; $toState<$numStates ; ++$toState)
  {
    my $line=<HMM>;
    my @fields=split/\s+/,$line;
    for(my $fromState=0 ; $fromState<$numStates ; ++$fromState)
      {
	$trans{$fromState}->{$toState}=$fields[$fromState];
      }
  }

# SKIP BLANK LINE
<HMM>;

# READ EMISSION PROBABILITIES
<HMM>;
my %emit;
foreach my $symbol (@alphabet)
  {
    my $line=<HMM>;
    my @fields=split/\s+/,$line;
    for(my $fromState=0 ; $fromState<$numStates ; ++$fromState)
      {
	$emit{$fromState}->{$symbol}=$fields[$fromState];
      }
  }

# OUTPUT TRANSITIONS
for(my $fromState=0 ; $fromState<$numStates ; ++$fromState)
  {
    my $line=<HMM>;
    my @fields=split/\s+/,$line;
    for(my $toState=0 ; $toState<$numStates ; ++$toState)
      {
	my $P=$trans{$fromState}->{$toState};
	next unless $P>0;
	print HMMS "$fromState -> $toState : $P\n";
      }
  }

# OUTPUT EMISSIONS
for(my $fromState=0 ; $fromState<$numStates ; ++$fromState)
  {
    print HMMS "state $fromState : ";
    foreach my $symbol (@alphabet)
      {
	my $P=0+$emit{$fromState}->{$symbol};
	if($fromState==0) {$P=0}
	print HMMS "$symbol=$P ";
      }
    print HMMS "\n";
  }

close(HMMS);
close(HMM);


