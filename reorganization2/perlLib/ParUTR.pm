package ParUTR;
use strict;
use ParSite;

######################################################################
#
# ParUTR.pm bmajoros@duke.edu 2/13/2012
#
# A UTR segment (i.e., untranslated portion of an exon)
#
# Attributes:
#   begin : 1-based inclusive coordinate of UTR begin
#   end : 1-based inclusive coordinate of UTR end
#   sites : array of ParSite
#   transcriptOffset : int
#
# Methods:
#   $utr=new ParUTR($begin,$end);
#   $utr->addSite($site);
#   $begin=$utr->getBegin();
#   $end=$utr->getEnd();
#   $len=$utr->getLength();
#   $numSites=$utr->getNumSites();
#   $site=$utr->getIthSite($i);
#   $bool=$utr->containsInterval($begin,$end);
#   $utr->setTranscriptOffset($offset);
#   $offset=$utr->getTranscriptOffset();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $utr=new ParUTR($begin,$end);
sub new
{
  my ($class,$begin,$end)=@_;
  my $self=
    {
     begin=>$begin,
     end=>$end,
     sites=>[]
    };
  bless $self,$class;
  return $self;
}
#---------------------------------------------------------------------
#   $utr->addSite($site);
sub addSite
  {
    my ($this,$site)=@_;
    push @{$this->{sites}},$site;
  }
#---------------------------------------------------------------------
#   $begin=$utr->getBegin();
sub getBegin
  {
    my ($this)=@_;
    return $this->{begin};
  }
#---------------------------------------------------------------------
#   $end=$utr->getEnd();
sub getEnd
  {
    my ($this)=@_;
    return $this->{end};
  }
#---------------------------------------------------------------------
#   $numSites=$utr->getNumSites();
sub getNumSites
  {
    my ($this)=@_;
    return 0+@{$this->{sites}};
  }
#---------------------------------------------------------------------
#   $site=$utr->getIthSite($i);
sub getIthSite
  {
    my ($this,$i)=@_;
    return $this->{sites}->[$i];
  }
#---------------------------------------------------------------------
#   $bool=$utr->containsInterval($begin,$end);
sub containsInterval
  {
    my ($this,$b,$e)=@_;
    return $b>=$this->{begin} && $e<=$this->{end};
  }
#---------------------------------------------------------------------
#   $utr->setTranscriptOffset($offset);
sub setTranscriptOffset
  {
    my ($this,$offset)=@_;
    $this->{transcriptOffset}=$offset;
  }
#---------------------------------------------------------------------
#   $offset=$utr->getTranscriptOffset();
sub getTranscriptOffset
  {
    my ($this)=@_;
    return $this->{transcriptOffset};
  }
#---------------------------------------------------------------------
#   $len=$utr->getLength();
sub getLength
  {
    my ($this)=@_;
    return $this->{end}-$this->{begin}+1;
  }
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

