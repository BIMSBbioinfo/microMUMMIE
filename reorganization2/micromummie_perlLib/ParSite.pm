package ParSite;
use strict;

######################################################################
#
# ParSite.pm bmajoros@duke.edu 2/13/2012
#
# Interaction site identified by PARalyzer (i.e., a "peak")
#
# Attributes:
#  begin : 1-based beginning coordinate (inclusive)
#  end : 1-based ending coordinate (inclusive)
#  score : float
#  libraryID : name protein
# Methods:
#   $parSite=new ParSite($begin,$end,$score,$libraryID);
#   $begin=$parSite->getBegin();
#   $end=$parSite->getEnd();
#   $length=$parSite->getLength();
#   $score=$parSite->getScore();
#   $lib=$parSite->getLibraryID();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $parSite=new ParSite($begin,$end,$score,$libraryID);
sub new
{
  my ($class,$begin,$end,$score,$libraryID)=@_;
  my $self=
    {
     begin=>$begin,
     end=>$end,
     score=>$score,
     libraryID=>$libraryID
    };
  bless $self,$class;
  return $self;
}
#---------------------------------------------------------------------
#   $begin=$parSite->getBegin();
sub getBegin
  {
    my ($this)=@_;
    return $this->{begin};
  }
#---------------------------------------------------------------------
#   $end=$parSite->getEnd();
sub getEnd
  {
    my ($this)=@_;
    return $this->{end};
  }
#---------------------------------------------------------------------
#   $score=$parSite->getScore();
sub getScore
  {
    my ($this)=@_;
    return $this->{score};
  }
#---------------------------------------------------------------------
#   $length=$parSite->getLength();
sub getLength
  {
    my ($this)=@_;
    return $this->{end}-$this->{begin};
  }
#---------------------------------------------------------------------
#   $lib=$parSite->getLibraryID();
sub getLibraryID
  {
    my ($this)=@_;
    return $this->{libraryID};
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

