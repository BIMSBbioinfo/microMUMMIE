#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;

my $name=ProgramName::get();
die "$name <in.gff> <groups.gff> <threshold>\n" unless @ARGV==3;
my ($infile,$groupsFile,$threshold)=@ARGV;

my $reader=new GffReader();
my $predictions=$reader->hashBySubstrate($infile);
my $groups=$reader->hashBySubstrate($groupsFile);
my @substrates=keys %$groups;
my $numSubstrates=@substrates;
for(my $i=0 ; $i<$numSubstrates ; ++$i) {
  my $substrate=$substrates[$i];
  my $preds=$predictions->{$substrate};
  my $grps=$groups->{$substrate};
  foreach my $pred (@$preds) {
    foreach my $group (@$grps){
      if($pred->overlapsOther($group)) {
	push @{$group->{sites}},$pred;
      }
    }
  }
  foreach my $group (@$grps) {
    next unless $group->{sites};
    my @sites=@{$group->{sites}};
    my $best;
    foreach my $site (@sites) {
      if(!$best || $site->{score}>$best->{score}) { $best=$site }
    }
    if($best->{score}<=$threshold) { print $best->toGff() }
  }
}






