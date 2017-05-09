#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <predictions.gff>\n" unless @ARGV==1;
my ($predFile)=@ARGV;

my $reader=new GffReader;
my $predHash=$reader->hashBySubstrate($predFile);
my $predictions=[];
my @substrates=keys %$predHash;
for(my $i=0 ; $i<@substrates ; ++$i) {
  my $substrate=$substrates[$i];
  my $preds=$predHash->{$substrate};
  @$preds=sort {
    $a->getBegin()<$b->getBegin() ? -1 :
      ($a->getBegin()>$b->getBegin() ? 1 :
        ($a->getEnd()<$b->getEnd() ? -1 :
          ($a->getEnd()>$b->getEnd() ? 1 : 0)
        )
      )
  } @$preds;
  my $n=@$preds;
  foreach my $prediction (@$preds) {
    my $add=join("",@{$prediction->{additionalFields}});
    if($add=~/_sh_/) {$prediction->{sh}=1}
  }

#goto SKIP;
  for(my $j=0 ; $j+1<$n ; ++$j) {
    my $first=$preds->[$j];
    my $last=$j;
    for(my $k=$j+1 ; $k<$n ; ++$k) {
      my $next=$preds->[$k];
      if($first->getBegin()!=$next->getBegin() ||
	 $first->getEnd()!=$next->getEnd()) { last }
      $last=$k;
    }
    next unless $last>$j;
    my ($numShuffled,$numNonShuffled);
    for(my $k=$j ; $k<=$last ; ++$k) 
      { if($preds->[$k]->{sh}) {++$numShuffled} else {++$numNonShuffled} }
    if($numShuffled>0 && $numNonShuffled>0) {
       my $numDel=$last-$j+1;
       splice(@$preds,$j,$numDel);
       $n-=$numDel;
       --$j;
    }
    else {
       my $numDel=$last-$j;
       splice(@$preds,$j+1,$numDel);
       $n-=$numDel;
       --$j;
    }
  }
#SKIP:
  for(my $j=0 ; $j<$n ; ++$j) {
    my $prediction=$preds->[$j];
    print $prediction->toGff();
  }
}
exit;

