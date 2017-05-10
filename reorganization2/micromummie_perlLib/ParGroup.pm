package ParGroup;
use strict;

######################################################################
#
# ParGroup.pm bmajoros@duke.edu 2/13/2012
#
# A "group" in PAR parlance (a scaffold of reads)
#
# Attributes:
#   begin : int
#   data : array of float
#   libraryID : string
#   signalType : string
# Methods:
#   $parGroup=new ParGroup($begin,$data,$libraryID,$signalType);
#   $begin=$group->getBegin();
#   $length=$group->getLength();
#   $data=$group->getData();
#   $libraryID=$group->getLibraryID();
#   $signalType=$group->getSignalType();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $parGroup=new ParGroup($begin,$data,$libraryID,$signalType);
sub new
{
  my ($class,$begin,$data,$libraryID,$signalType)=@_;
  my $self=
    {
     begin=>$begin,
     data=>$data,
     libraryID=>$libraryID,
     signalType=>$signalType
    };
  bless $self,$class;
  return $self;
}
#---------------------------------------------------------------------
#   $begin=$group->getBegin();
sub getBegin
  {
    my ($this)=@_;
    return $this->{begin};
  }
#---------------------------------------------------------------------
#   $length=$group->getLength();
sub getLength
  {
    my ($this)=@_;
    return 0+@{$this->{data}};
  }
#---------------------------------------------------------------------
#   $data=$group->getData();
sub getData
  {
    my ($this)=@_;
    return $this->{data};
  }
#---------------------------------------------------------------------
#   $libraryID=$group->getLibraryID();
sub getLibraryID
  {
    my ($this)=@_;
    return $this->{libraryID};
  }
#---------------------------------------------------------------------
#   $signalType=$group->getSignalType();
sub getSignalType
  {
    my ($this)=@_;
    return $this->{signalType};
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

