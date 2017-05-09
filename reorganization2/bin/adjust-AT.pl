#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <background-AT> <group-AT>\n" unless @ARGV==2;
my ($bgAT,$flankAT)=@ARGV;

if($bgAT>=1 || $flankAT>=1) { die "probabilities must be between 0 and 1\n" }

process("bg.hmm",$bgAT);
process("flank-trained.hmm",$flankAT);

sub process
  {
    my ($file,$AT)=@_;
    my $halfAT=$AT/2;
    my $halfGC=(1-$AT)/2;
    system("bin/hmm-edit $file NMER 1 0 A $halfAT NMER 1 0 T $halfAT NMER 1 0 G $halfGC NMER 1 0 C $halfGC");
  }


