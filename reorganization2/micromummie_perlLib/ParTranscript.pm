package ParTranscript;
use strict;
use ParUTR;
use ParGroup;

######################################################################
#
# ParTranscript.pm bmajoros@duke.edu 2/13/2012
#
# Attributes:
#   transcriptID : string
#   UTRs : array of ParUTR
# Methods:
#   $parTranscript=new ParTranscript($transcriptID);
#   $id=$parTranscript->getID();
#   $parTranscript->addUTR($utr);
#   $n=$parTranscript->numUTRs();
#   $utr=$parTranscript->getIthUTR($i);
#   $parTranscript->sortUTRs($strand);
#   $bool=$transcript->isLongerThan($L);
#   $bool=$transcript->overlapsGroup($parGroup);
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $parTranscript=new ParTranscript($transcriptID);
sub new
{
  my ($class,$transcriptID)=@_;
  my $self=
    {
     transcriptID=>$transcriptID,
     UTRs=>[]
    };
  bless $self,$class;
  return $self;
}
#---------------------------------------------------------------------
#   $parTranscript->addUTR($utr);
sub addUTR
  {
    my ($this,$utr)=@_;
    push @{$this->{UTRs}},$utr;
  }
#---------------------------------------------------------------------
#   $n=$parTranscript->numUTRs();
sub numUTRs
  {
    my ($this)=@_;
    return 0+@{$this->{UTRs}};
  }
#---------------------------------------------------------------------
#   $utr=$parTranscript->getIthUTR($i);
sub getIthUTR
  {
    my ($this,$i)=@_;
    return $this->{UTRs}->[$i];
  }
#---------------------------------------------------------------------
#   $id=$parTranscript->getID();
sub getID
  {
    my ($this)=@_;
    return $this->{transcriptID};
  }
#---------------------------------------------------------------------
#   $parTranscript->sortUTRs($strand);
sub sortUTRs
  {
    my ($this,$strand)=@_;
    my $UTRs=$this->{UTRs};
    if($strand eq "+")
      { @$UTRs=sort {$a->getBegin() <=> $b->getBegin()} @$UTRs }
    elsif($strand eq "-")
      { @$UTRs=sort {$b->getBegin() <=> $a->getBegin()} @$UTRs }
    else { die "bad strand: $strand" }
    my $n=@$UTRs;
    my $sum=0;
    for(my $i=0 ; $i<$n ; ++$i) {
      my $utr=$UTRs->[$i];
      $utr->setTranscriptOffset($sum);
      $sum+=$utr->getLength();
    }
  }
#---------------------------------------------------------------------
#   $bool=$transcript->isLongerThan($L);
sub isLongerThan
  {
    my ($this,$L)=@_;
    my $UTRs=$this->{UTRs};
    my $n=@$UTRs;
    if($n<1) {die}
    my $min=$UTRs->[0]->getBegin();
    my $max=$UTRs->[0]->getEnd();
    for(my $i=1 ; $i<$n ; ++$i) {
      my $utr=$UTRs->[$i];
      my $begin=$utr->getBegin();
      my $end=$utr->getEnd();
      if($begin<$min) {$min=$begin}
      if($end>$max) {$max=$end}
    }
    return $max-$min>$L;
  }
#---------------------------------------------------------------------
#   $bool=$transcript->overlapsGroup($parGroup);
sub overlapsGroup
  {
    my ($this,$group)=@_;
    my $groupBegin=$group->getBegin();
    my $groupEnd=$groupBegin+$group->getLength()-1;
    my $UTRs=$this->{UTRs};
    foreach my $utr (@$UTRs) {
      if($utr->containsInterval($groupBegin,$groupEnd)) { return 1 }
    }
    return 0;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

