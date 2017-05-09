#!/usr/bin/perl
use strict;
use TempFilename;
use File::Basename;
use File::Spec;
use File::chdir;

my $VERBOSE=1;

my $WANT_BULGE=0;
my $WANT_CONSERVATION=0;
my $TARGETSCAN_VARIANCE=3.0;
my $TARGETSCAN_BG_MEAN=0.0;
my $TARGETSCAN_FG_MEAN=3.0;
my $dir_path = dirname(File::Spec->rel2abs(__FILE__));

die "microMUMMIE.pl <mature-miRNAs.txt> <genome.2bit> <paralyzer-output-dir> <library-name> <out.gff> <posterior-decoding:0/1> <UTRs.txt> <Path of directory where you want to write all analysis in the end>\n"
  unless @ARGV==8;
my ($mature,$twoBitFile,$dataDir,$libraryName,$outfile,$wantPost,$UTR_FILE,$dir)=@ARGV;
my $DASH_G=$wantPost ? "-I" : "-g";

System("date");

`cp $UTR_FILE site.hmm bg.hmm peak.hmm flank-trained.hmm metamodel.txt submodels.txt cons.schema nocons.schema make-tgf.pl hmm-edit random-HMM fasta-to-fastb.pl baum-welch model-combiner $dir`;

my $groupsFile="$dataDir/$libraryName.groups.csv";

$CWD=$dir;
print "
# -------------------------------------------
# EXTRACTING SEEDS FROM MATURE MICRO-RNA LIST
# -------------------------------------------
";

System("$dir_path/get-seeds.pl $mature seeds.fasta");

print "
# -------------------------------------------
# BUILDING THE SITE SUBMODEL
# -------------------------------------------
";

System("$dir_path/make-42state-site.pl seeds.fasta .33 .33 .33 4 cons.schema 7 bg.hmm site.hmm '' N N");

print "
# -------------------------------------------
# COMBINING SUBMODELS INTO FULL MODEL
# -------------------------------------------
";

System("model-combiner metamodel.txt submodels.txt PARCLIP.hmm");
if($WANT_CONSERVATION)
  {
  System("hmm-edit PARCLIP.hmm MEAN all -- 1 $TARGETSCAN_BG_MEAN");
  System("hmm-edit PARCLIP.hmm VAR  all -- 1 $TARGETSCAN_VARIANCE");
  System("hmm-edit PARCLIP.hmm MEAN 3   -- 1 $TARGETSCAN_FG_MEAN");
  }
else {
    System("hmm-edit PARCLIP.hmm DTRK targetscan");
}


print "
# -------------------------------------------
# PREPARING INPUT FILES
# -------------------------------------------
";

System("rm -rf $dir/chunks") if -e "chunks";
System("rm -rf $dir/chunks2") if -e "chunks2";

System("mkdir $dir/chunks $dir/chunks2");

System("$dir_path/par2fastb.pl $twoBitFile $dataDir $dir/chunks $UTR_FILE $libraryName");

System("$dir_path/assemble-transcripts.pl chunks chunks2 ; rm -r $dir/chunks ; mv $dir/chunks2 $dir/chunks");

print "
# -------------------------------------------
# RUN THE MODEL
# -------------------------------------------
";

System("rm predictions-var*.gff");
my @vars=(0.5,   0.25,  0.2,  0.15,  0.1,  0.01);
my @sens=(0.12,  0.17,  0.2,  0.27,  0.42, 0.62);
my @SNR= (15.7, 12.04,  9.95, 7.07,  5.09, 2.24);
my $N=@vars;
for(my $i=0 ; $i<$N ; ++$i) {
  my $var=$vars[$i];  my $sens=$sens[$i];  my $snr=$SNR[$i];
  print "Running at variance $var\n";
  System("$dir_path/hmm-edit PARCLIP.hmm VAR all 0 $var");
  System("$dir_path/parse $DASH_G 5-45 -p -d PARCLIP.hmm $dir/chunks > chunk-preds.gff");
  System("$dir_path/identify-miRNAs.pl chunk-preds.gff seeds.fasta > identify.tmp");
  System("$dir_path/combine-miRNA-predictions.pl identify.tmp > chunk-preds.gff");
  System("$dir_path/get-chrom-coords_try.pl chunk-preds.gff > predictions-var$var.gff");
  addScores("predictions-var$var.gff",$sens,$snr);
  System("$dir_path/match_coordinates1.pl predictions-var$var.gff chunks predictions-var$var-genomic.gff");
  `cut -f 1,4,5,6,7,8,9 predictions-var$var-genomic.gff > test.gff`;
  `cat $groupsFile | perl -wp -i.backup -ne 's/,/\t/g' | while read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 col16 col17 col18 col19 col20 col21 col22; do echo \$col1 \$col3 \$col4 \$col8 \$col10 \$col2 \$col5; done | perl -wp -i.backup -ne 's/ /\t/g' | perl -ne'\$.==1?next:print' | intersectBed -wao -a test.gff -b stdin | perl -e 'while(<>){ \@line = split(/\t/, \$_); if( \$line[14]=~ 0 ) { next;} else  { print "\$_"; } }' | while read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15; do echo \$col1 \$col14 \$col11 \$col2 \$col3 \$col4 \$col13 \$col7 \$col9 \$col10; done | perl -i -ple 'print q{Chromosome Alinged_to GroupID Prediction_start Prediction_end Score Strand Info Group_start Group_end} if \$. == 1; close ARGV if eof' | perl -wp -i.backup -ne 's/ /\t/g' > predictions-var$var-map.gff`;
`cut -f 2 predictions-var$var-map.gff | sort | uniq -c > alignment_stats_var$var`;
  if($WANT_BULGE) {
    print "Running bulge models\n";
    my $type="";
    for(my $bulgeType=1 ; $bulgeType<=3 ; ++$bulgeType) {
      $type.="I";
      System("$dir_path/bulge.pl chunks seeds.fasta $bulgeType $var");
      System("$dir_path/get-chrom-coords_try.pl bulge-predictions-var$var.gff > chunk-preds.gff");
      System("$dir_path/mv chunk-preds.gff bulge$type-predictions-var$var.gff");
      addScores("bulge$type-predictions-var$var.gff",$sens,$snr);
      System("$dir_path/match_coordinates1.pl bulge$type-predictions-var$var.gff chunks bulge$type-predictions-var$var-genomic.gff");
      `cut -f 1,4,5,6,7,8,9 bulge$type-predictions-var$var-genomic.gff > test.gff`;
      `cat $groupsFile | perl -wp -i.backup -ne 's/,/\t/g' | while read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 col16 col17 col18 col19 col20 col21 col22; do echo \$col1 \$col3 \$col4 \$col8 \$col10 \$col2 \$col5; done | perl -wp -i.backup -ne 's/ /\t/g' | perl -ne'\$.==1?next:print' | intersectBed -wao -a test.gff -b stdin | perl -e 'while(<>){ \@line = split(/\t/, \$_); if( \$line[14]=~ 0 ) { next;} else  { print "\$_"; } }' | while read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15; do echo \$col1 \$col14 \$col11 \$col2 \$col3 \$col4 \$col13 \$col7 \$col9 \$col10; done | perl -i -ple 'print q{Chromosome Alinged_to GroupID Prediction_start Prediction_end Score Strand Info Group_start Group_end} if \$. == 1; close ARGV if eof' | perl -wp -i.backup -ne 's/ /\t/g' > bulge$type-predictions-var$var-map.gff`;
`cut -f 2 bulge$type-predictions-var$var-map.gff | sort | uniq -c > alignment_stats_bulge_var$var`;
    }
  }
}

`ls predictions-var*.gff > filename`;
`cat filename | grep -v genomic | grep -v map > filename1`;
`cat filename | grep genomic > filename2`;
`cat filename | grep map > filename3`;
System("for fn in `cat filename1`; do cat \$fn; done > $outfile");
System("for fn in `cat filename2`; do cat \$fn | perl -ni.bak -e'print unless m/^Chromosome/' ; done > $outfile-genomic.gff");
System("for fn in `cat filename3`; do cat \$fn | perl -ni.bak -e'print unless m/^Chromosome/' ; done > $outfile-map.gff");
`cut -f 2 $outfile-map.gff | sort | uniq -c > alignment_stats_$outfile-genomic`;


print "
# -------------------------------------------
# DONE.  OUTPUT IS IN: $outfile
# -------------------------------------------
";

System("date");

`rm out.gff test.gff chunk-preds.gff site.hmm bg.hmm peak.hmm flank-trained.hmm metamodel.txt submodels.txt cons.schema nocons.schema make-tgf.pl hmm-edit random-HMM fasta-to-fastb.pl baum-welch model-combiner $UTR_FILE filename filename1 filename2 filename3`;

sub addScores
  {
    my ($filename,$sens,$snr)=@_;
    my $tempName=TempFilename::generate();
    open(IN,$filename) || die "can't open file: $filename\n";
    open(OUT,">$tempName") || die "can't write to file: $tempName\n";
    while(<IN>) {
      chomp;
      print OUT "${_}sens=$sens;SNR=$snr;\n";
    }
    close(OUT);
    close(IN);
    System("mv $tempName $filename");
  }


sub System
{
    my ($cmd)=@_;
    if($VERBOSE) { print "$cmd\n" }
    system($cmd);
}

