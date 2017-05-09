#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use TempFilename;

#================================================================
#   miRNA index:    87654321
#          6mer3-8: XXXXXX--
#          6mer2-7: -XXXXXX-
#          7mer-m8: XXXXXXX-
#          7mer-m1: -XXXXXXX
#             8mer: XXXXXXXX
#================================================================

my $name=ProgramName::get();
die "$name <8mers.fasta> <6mer-prob> <7mer-prob> <8mer-prob> <num-mix-comp> <schema> <order> <bg.hmm> <out.hmm> <weights.txt> <TS_seed_types:Y/N> <uniform:Y/N>\n" unless @ARGV==12;
my ($seedsFile,$sixProb,$sevenProb,$eightProb,$mix,$schema,$order,$bgHMM,
    $outfile,$weightFile,$TS_seed_types,$uniformSeeds)=@ARGV;
$TS_seed_types="\U$TS_seed_types";
$TS_seed_types=($TS_seed_types eq "Y");

# Misc initialization
my $SUBDIR=TempFilename::generate();
system("mkdir $SUBDIR");
die unless -e $SUBDIR;
system("cp $bgHMM $SUBDIR/bg.hmm");
writeMetamodel();
my $tgfFile="$SUBDIR/tgf.tgf";
system("make-tgf.pl $schema full $tgfFile");
makeAfiles("$SUBDIR/A.hmm",$schema);

my $DASH_W;
if($weightFile ne "") { $DASH_W="-w $weightFile" }

# Load 8mers
my @seeds;
my $reader=new FastaReader($seedsFile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless $def;
  push @seeds,[$def,$seq];
}

# Make canonical models for each length
makeSubmodel(0,5,"6mer3-8");#          6mer3-8: XXXXXX--
makeSubmodel(1,6,"6mer2-7");#          6mer2-7: -XXXXXX-
makeSubmodel(0,6,"7mer-m8");#          7mer-m8: XXXXXXX-
makeSubmodel(1,7,"7mer-m1");#          7mer-m1: -XXXXXXX
makeSubmodel(0,7,"8mer");   #             8mer: XXXXXXXX

# Combine into complete model
my $pwd=`pwd`;
chomp $pwd;
system("cd $SUBDIR ; model-combiner metamodel.txt submodels.txt $pwd/$outfile");

# Set mixture parameters
system("hmm-edit $outfile MIX ALL 0 0 MIX ALL 1 0 MIX ALL 2 0 MIX ALL 3 1");
system("hmm-edit $outfile MEAN 0 0 1e-7   MEAN 0 1 0.531313");
system("hmm-edit $outfile MEAN 1 0 0.5    MEAN 1 1 0.531313");
system("hmm-edit $outfile MEAN 2 0 1.0    MEAN 2 1 0.531313");
system("hmm-edit $outfile MEAN 3 0 0.5    MEAN 3 1 0.993844");
system("hmm-edit $outfile COV ALL ALL ALL 0");
system("hmm-edit $outfile VAR  0 0 0.01   VAR  0 1 1.10313");
system("hmm-edit $outfile VAR  1 0 0.01   VAR  1 1 1.10313");
system("hmm-edit $outfile VAR  2 0 0.01   VAR  2 1 1.10313");
system("hmm-edit $outfile VAR  3 0 0.01   VAR  3 1 1.15067");
if($TS_seed_types) {
  print STDERR "using TS site types\n";
  system("hmm-edit $outfile TRANS 0 34 0.0 TRANS 0 26 0.0");
  system("hmm-edit $outfile TRANS 0 18 0.333 TRANS 0 10 0.333 TRANS 0 1 0.333");
  system("hmm-edit $outfile TRANS 7 9 1.0 TRANS 16 9  1.0");
  system("hmm-edit $outfile TRANS 7 8 0.0 TRANS 16 17 0.0");
}
elsif($uniformSeeds) {
  my $p=1/7.0;
  my $p2=2/7.0;
  system("hmm-edit $outfile TRANS 0 1 $p2 TRANS 0 10 $p2 TRANS 0 18 $p TRANS 0 26 $p TRANS 0 34 $p");
  system("hmm-edit $outfile TRANS 7 9 0.5 TRANS 16 9 0.5"); # 1A forms
  system("hmm-edit $outfile TRANS 7 8 0.5 TRANS 16 17 0.5"); # 1A forms
}
else {
  system("hmm-edit $outfile TRANS 7 9 0.5 TRANS 16 9 0.5"); # 1A forms
  system("hmm-edit $outfile TRANS 7 8 0.5 TRANS 16 17 0.5"); # 1A forms
}
# Clean up
system("rm -r $SUBDIR");


######################################################################3
sub makeSubmodel {
  my ($begin,$end,$type)=@_;
  my $len=$end-$begin+1;
  my $hmmFile="$SUBDIR/$type.hmm";
  my $trainFile="$SUBDIR/$type.seeds";
  open(OUT,">$trainFile") || die;
  foreach my $seed (@seeds) {
    my ($def,$seq)=@$seed;
    $def=~/>(\S+)/ || die "$def $seed $begin $end $type";
    my $id=$1;
    my $subseed=substr($seq,$begin,$len);
    print OUT ">$id\n$subseed\n";
  }
  close(OUT);
  my $nextSeed=int(rand(1000000));
  my $numStates=$len+1;
  my $order=$len-1;
  my $templateHMM="$hmmFile.template";
  system("random-HMM -c V -u -s $nextSeed $numStates 1 $mix $schema $order $templateHMM");
  my $tieProfile="$SUBDIR/$type.tie";
  makeTieProfile($tieProfile,$len);
  my $nmersDir="$SUBDIR/$type.nmers";
  system("mkdir $nmersDir");
  system("fasta-to-fastb.pl $trainFile $schema $nmersDir");
  my $nextSeed=int(rand(1000000));
  system("baum-welch $templateHMM $tgfFile $nmersDir 2 $hmmFile -R -s $nextSeed -u -t $tieProfile -n 0.00001 $DASH_W");
  system("rm $trainFile $tieProfile $templateHMM");
  system("rm -r $nmersDir");
}
######################################################################3
sub writeMetamodel {
  my $metamodel="$SUBDIR/metamodel.txt";
  my $submodels="$SUBDIR/submodels.txt";
  my $halfSeven=$sevenProb/2;
  my $halfSix=$sixProb/2;
  open(OUT,">$submodels") || die;
  print OUT 
"1 = 8mer.hmm
2 = A.hmm
3 = bg.hmm
4 = 7mer-m1.hmm
5 = 7mer-m8.hmm
6 = bg.hmm
7 = bg.hmm
8 = 6mer2-7.hmm
9 = bg.hmm
10 = 6mer3-8.hmm
11 = bg.hmm
12 = bg.hmm
";
  close(OUT);
  open(OUT,">$metamodel") || die;
  print OUT 
"0 -> 1 : $eightProb
0 -> 3 : $halfSeven
0 -> 5 : $halfSeven
0 -> 7 : $halfSix
0 -> 10 : $halfSix
1 -> 0 : 1
2 -> 0 : 1
3 -> 4 : 1
4 -> 0 : 1
5 -> 6 : 1
6 -> 0 : 1
7 -> 8 : 1
8 -> 9 : 1
9 -> 0 : 1
10 -> 11 : 1
11 -> 12 : 1
12 -> 0 : 1
";
  close(OUT);
}
######################################################################3
sub makeAfiles {
  my ($outHMM,$schema)=@_;
  my $fasta=TempFilename::generate();
  my $dir=TempFilename::generate();
  system("mkdir $dir");
  open(OUT,">$fasta") || die;
  print OUT ">1\nA\n>2\nA\n>3\nA\n";
  close(OUT);
  system("fasta-to-fastb.pl $fasta $schema $dir");
  my $inHMM=TempFilename::generate();
  system("random-HMM -c 0 2 1 $mix $schema $order $inHMM");
  system("hmm-edit $inHMM TRANS 1 1 0");
  my $tieProfile=TempFilename::generate();
  makeTieProfile($tieProfile,1);
  my $nextSeed=int(rand(1000000));
  system("baum-welch $inHMM $tgfFile $dir 2 $outHMM -R -s $nextSeed -u -t $tieProfile -n 0.00001");
  system("rm $fasta $inHMM $tieProfile ; rm -r $dir");
}
######################################################################3
sub makeTieProfile {
  my ($filename,$L)=@_;
  open(OUT,">$filename") || die $filename;
  print OUT "fix means\nfix covariance_matrix\nfix transitions\n";
  print OUT "fix weights in states ";
  for(my $q=1 ; $q<$L ; ++$q) { print OUT "$q," }
  print OUT "$L\n";
  close(OUT);
}
######################################################################3
######################################################################3
######################################################################3
######################################################################3
######################################################################3



