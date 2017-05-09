#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;
$|=1;
my $CWD = `pwd`;        chomp $CWD;  $CWD =~ s/ //gi;
my $SCHEMA="cons.schema";
my $name=ProgramName::get();
die "$name <chunks-dir> <seeds.fasta> <rundir> <bulge-type> <var> \n" 
  unless @ARGV==5;
my ($TEST_DIR,$UNSHUFFLED,$RUNDIR,$BULGE_TYPE,$VARIANCE)=@ARGV;
system("mkdir -p $RUNDIR");
my $DASH_G="-g ";  #already has -p down below  #"-g"; # "-g" or "-I";
my $TRAINING_SEEDS="seeds.fasta";
my $METAMODEL="bulge/bulge-metamodel.txt";#"shuffled-metamodel.txt"; # N3
my $SUBMODELS="bulge/bulge-submodels$BULGE_TYPE.txt";#"shuffled-submodels.txt"; # N3
my $STATES="5-11"; #"5-45"; # N3 w/42-state site
my $PEAK_THRESHOLD=0.999;

# Load the unshuffled seeds
my @unshuffled;
my $reader=new FastaReader($UNSHUFFLED);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  $def=~/>(\S+)/ || die;
  my $id=$1;
  push @unshuffled,[$id,$seq];
}
$reader->close();

my $X="$RUNDIR/$BULGE_TYPE";

my @sample;
foreach my $pair (@unshuffled) {  push @sample,$pair; }

open(OUT,">$X.$TRAINING_SEEDS") || die;
foreach my $pair (@sample) {
  my ($def,$seq)=@$pair;
  print OUT ">$def\n$seq\n";
}
close(OUT);

system("
perl/nicky/make-bulge-site.pl $X.$TRAINING_SEEDS 4 cons.schema 6 $X.site.hmm 
cp $SUBMODELS $X.submodels
perl/sub.pl $X.submodels bulge-site$BULGE_TYPE.hmm $X.site.hmm
bin/model-combiner $METAMODEL $X.submodels $X.hmm
rm $X.submodels
bin/hmm-edit $X.hmm DTRK phastcons
bin/hmm-edit $X.hmm VAR all 0 $VARIANCE
bin/parse $X.hmm $TEST_DIR -d -p $DASH_G $STATES > $X.predictions.gff.tmp
bulge/identify-miRNAs-bulge.pl $X.predictions.gff.tmp $X.$TRAINING_SEEDS > $X.predictions.gff.tmp2
rm $X.predictions.gff.tmp
perl/combine-miRNA-predictions.pl $X.predictions.gff.tmp2 > $X.predictions.gff
rm $X.predictions.gff.tmp2
#cp $X.predictions.gff tmp.gff
bulge/identify-miRNAs-bulge.pl $X.predictions.gff $X.$TRAINING_SEEDS > bulge-predictions-var$VARIANCE.gff
");


