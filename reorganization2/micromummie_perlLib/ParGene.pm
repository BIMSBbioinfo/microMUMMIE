package ParGene;
use strict;
use ParTranscript;
use ParGroup;
######################################################################
#
# ParGene.pm bmajoros@duke.edu 2/13/2012
#
# Attributes:
#   chr : string
#   strand : '+' or '-'
#   transcripts : hash mapping transcriptID to ParTranscript
#   groups : array of ParGroup
# Methods:
#   $parGene=new ParGene($chr,$strand);
#   $chr=$parGene->getChr();
#   $strand=$parGene->getStrand();
#   $transcript=$parGene->findOrCreateTranscript($transcriptID);
#   $transcript=$parGene->find($transcriptID);
#   $transcriptIDs=$parGene->getTranscriptIDs();
#   $parGene->addGroup($parGroup);
#   $groups=$gene->getGroups();
#   $bool=$gene->isLongerThan($L);
#   $gene->removeTranscript($transcriptID);
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $parGene=new ParGene($chr,$strand);
sub new
{
  my ($class,$chr,$strand)=@_;
  my $self=
    {
     chr=>$chr,
     strand=>$strand,
     transcripts=>{}
    };
  bless $self,$class;
  return $self;
}
#---------------------------------------------------------------------
#   $chr=$parGene->getChr();
sub getChr
  {
    my ($this)=@_;
    return $this->{chr};
  }
#---------------------------------------------------------------------
#   $strand=$parGene->getStrand();
sub getStrand
  {
    my ($this)=@_;
    return $this->{strand};
  }
#---------------------------------------------------------------------
#   $transcript=$parGene->findOrCreateTranscript($transcriptID);
sub findOrCreateTranscript
  {
    my ($this,$transcriptID)=@_;
    my $hash=$this->{transcripts};
    if(!defined($hash->{$transcriptID}))
      { $hash->{$transcriptID}=new ParTranscript($transcriptID) }
    return $hash->{$transcriptID};
  }
#---------------------------------------------------------------------
#   $transcript=$parGene->find($transcriptID);
sub find
  {
    my ($this,$transcriptID)=@_;
    my $hash=$this->{transcripts};
    #if(!defined($hash->{$transcriptID})) {die "$transcriptID not found"}
    return $hash->{$transcriptID};
  }
#---------------------------------------------------------------------
#   $transcriptIDs=$parGene->getTranscriptIDs();
sub getTranscriptIDs
  {
    my ($this)=@_;
    my $array=[];
    @$array=keys %{$this->{transcripts}};
    return $array;
  }
#---------------------------------------------------------------------
#   $parGene->addGroup($parGroup);
sub addGroup
  {
    my ($this,$group)=@_;
    push @{$this->{groups}},$group;
  }
#---------------------------------------------------------------------
#   $groups=$gene->getGroups();
sub getGroups
  {
    my ($this)=@_;
    return $this->{groups};
  }
#---------------------------------------------------------------------
#   $bool=$gene->isLongerThan($L);
sub isLongerThan
  {
    my ($this,$L)=@_;
    my $transcriptIDs=$this->getTranscriptIDs();
    foreach my $transcriptID (@$transcriptIDs) {
      my $transcript=$this->{transcripts}->{$transcriptID};
      if($transcript->isLongerThan($L)) { return 1 }
    }
    return 0;
  }
#---------------------------------------------------------------------
#   $gene->removeTranscript($transcriptID);
sub removeTranscript
  {
    my ($this,$id)=@_;
    my $transcripts=$this->{transcripts};
    delete $transcripts->{$id};
  }
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

