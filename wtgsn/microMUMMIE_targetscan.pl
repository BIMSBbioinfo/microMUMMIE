#!/usr/bin/perl
use strict;
use TempFilename;
use File::Basename;
use File::Spec;
use File::chdir;

my $VERBOSE=0;
my $WANT_BULGE=0;
my $WANT_CONSERVATION=1;
my $TARGETSCAN_VARIANCE=3.0;
my $TARGETSCAN_BG_MEAN=0.0;
my $TARGETSCAN_FG_MEAN=3.0;
my $TARGETSCAN_SCORETYPE="branchLen";
my $TOTAL_ALIGNMENT_FILE="mchr_orderedtc_scUTRs.maf";
my $MAPPED_ID_FILE="mappedspid23w.txt";
my $dir_path = dirname(File::Spec->rel2abs(__FILE__));
my $mummie_path="/data/ohler/Software/microMUMMIE_new/wtgsn";



die "microMUMMIE_targetscan.pl <mature-miRNAs.txt> <mature.fa> <genome.2bit> <paralyzer-output-dir> <library-name> <out.gff> <posterior-decoding:0/1> <UTRs.txt> <Maf file> <Path of directory where you want to write all analysis in the end>\n"
  unless @ARGV==10;
my ($mature,$maturefa,$twoBitFile,$dataDir,$libraryName,$outfile,$wantPost,$UTR_FILE,$maf,$dir)=@ARGV;
my $DASH_G=$wantPost ? "-I" : "-g";

my $groupsFile="$dataDir/$libraryName.groups.csv";
System("date");

print "
# -------------------------------------------
# PREPARING FILES
# -------------------------------------------
";

`cp -r $mummie_path/pct_trees_parameters $dir`;
`cp pct_trees_parameters $UTR_FILE $mummie_path/*.hmm $mummie_path/*.pl $mummie_path/*.pm $mummie_path/metamodel.txt $mummie_path/submodels.txt $mummie_path/cons.schema $mummie_path/nocons.schema $mummie_path/baum-welch $mummie_path/classify $mummie_path/extract-chain $mummie_path/extract-motifs $mummie_path/get-likelihood $mummie_path/hmm-edit $mummie_path/hmm-extract-state $mummie_path/install-chain $mummie_path/kmeans $mummie_path/model-combiner $mummie_path/parse $mummie_path/random-HMM $mummie_path/sample $mummie_path/sub $mummie_path/xgraph $dir`;



System("rm -rf $dir/chunks") if -e "$dir/chunks";
System("rm -rf $dir/chunks-targetscan") if -e "$dir/chunks-targetscan";
System("rm -rf $dir/chunks-assembled") if -e "$dir/chunks-assembled";

System("mkdir $dir/chunks $dir/chunks-targetscan $dir/chunks-assembled");


$CWD=$dir;

System("$dir_path/par2fastb.pl $twoBitFile $dataDir $dir/chunks $UTR_FILE $libraryName");

System("$dir_path/assemble-transcripts.pl $dir/chunks $dir/chunks-assembled");

$CWD = "$dir/chunks-assembled";

`grep ">DNA" * > $dir/file_defid.txt`;

$CWD = $dir;

System("$dir_path/get_id_from_chunks.pl");

$CWD=$dir_path;

System("perl $dir_path/fastamiR_to_targetscan_input.pl $maturefa $dir/miRNA_txt");

print "
# -------------------------------------------
# EXTRACTING ALIGNMENT FOR TARGETSCAN
# -------------------------------------------
 ";

System("$dir_path/parseandconcatalmfortgscn.pl $maf $MAPPED_ID_FILE $dir/out_enst.txt $dir/out.maf" );


print "
# -------------------------------------------
# RUNNING TARGETSCAN NOW
# -------------------------------------------
";

$CWD=$dir;

System("perl $dir_path/targetscan_60.pl miRNA_txt out.maf targetscan_60_output.txt");

$CWD=$dir_path;

System("$dir_path/targetscan_60_BL_bins.pl $dir/out.maf > $dir/UTRs_median_BLs_bins.my_output.txt");

$CWD=$dir;

if (-e "targetscan_60_output.sort.txt") {
System("rm targetscan_60_output.sort.txt") 
}

$CWD=$dir_path;

System("perl $dir_path/targetscan_60_BL_PCT.pl $dir/miRNA_txt $dir/targetscan_60_output.txt $dir/UTRs_median_BLs_bins.my_output.txt > $dir/targetscan_60_output.BL_PCT.my_output.txt");

System("perl $dir_path/targetscan_60_context_scores.pl $dir/miRNA_txt_cs $dir/out.maf $dir/targetscan_60_output.txt $dir/targetscan_60_context_scores_output.txt");

$CWD=$dir;

print "
# -------------------------------------------
# # # ADDING TARGETSCAN TRACK TO FASTB FILES
# # # -------------------------------------------
 ";

System("$dir_path/formattgscnoutputforbill.pl targetscan_60_context_scores_output.txt targetscan_60_output.BL_PCT.my_output.txt tgscn_combined_out.txt 9606");

System("$dir_path/get-transcript-ids.pl chunks-assembled > transcript-ids.txt");

System("$dir_path/get-targetscan.pl tgscn_combined_out.txt $maturefa Y transcript-ids.txt 1 > targetscan.gff");

System("$dir_path/add-targetscan-tracks.pl targetscan.gff chunks-assembled chunks-targetscan $TARGETSCAN_SCORETYPE");

System("$dir_path/assemble-transcripts.pl chunks-targetscan chunks-assembled");




print "
# -------------------------------------------
# EXTRACTING SEEDS FROM MATURE MICRO-RNA LIST
# -------------------------------------------
";

System("perl get-seeds.pl $mature seeds.fasta");


print "
# -------------------------------------------
# BUILDING THE SITE SUBMODEL
# -------------------------------------------
";

System("make-42state-site.pl seeds.fasta .33 .33 .33 4 cons.schema 7 bg.hmm site.hmm '' N N");

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
  System("hmm-edit PARCLIP.hmm VAR all 0 $var");
  System("parse $DASH_G 5-45 -p -d PARCLIP.hmm chunks-targetscan > chunk-preds.gff");
  System("identify-miRNAs.pl chunk-preds.gff seeds.fasta > identify.tmp");
  System("combine-miRNA-predictions.pl identify.tmp > chunk-preds.gff");
  System("get-chrom-coords_try.pl chunk-preds.gff > predictions-var$var.gff");
  addScores("predictions-var$var.gff",$sens,$snr);
  System("$dir_path/match_coordinates1.pl predictions-var$var.gff chunks-assembled predictions-var$var-genomic.gff");
  `cut -f 1,4,5,6,7,8,9 predictions-var$var-genomic.gff > test.gff`;
  `cat $groupsFile | perl -wp -i.backup -ne 's/,/\t/g' | while read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15 col16 col17 col18 col19 col20 col21 col22; do echo \$col1 \$col3 \$col4 \$col8 \$col10 \$col2 \$col5; done | perl -wp -i.backup -ne 's/ /\t/g' | perl -ne'\$.==1?next:print' | intersectBed -wao -a test.gff -b stdin | perl -e 'while(<>){ \@line = split(/\t/, \$_); if( \$line[14]=~ 0 ) { next;} else  { print "\$_"; } }' | while read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14 col15; do echo \$col1 \$col14 \$col11 \$col2 \$col3 \$col4 \$col13 \$col7 \$col9 \$col10; done | perl -i -ple 'print q{Chromosome Alinged_to GroupID Prediction_start Prediction_end Score Strand Info Group_start Group_end} if \$. == 1; close ARGV if eof' | perl -wp -i.backup -ne 's/ /\t/g' > predictions-var$var-map.gff`;
`cut -f 2 predictions-var$var-map.gff | sort | uniq -c > alignment_stats_var$var`;
  if($WANT_BULGE) {
    print "Running bulge models\n";
    my $type="";
    for(my $bulgeType=1 ; $bulgeType<=3 ; ++$bulgeType) {
      $type.="I";
      System("bulge.pl chunks-targetscan seeds.fasta $bulgeType $var");
      System("get-chrom-coords_try.pl bulge-predictions-var$var.gff > chunk-preds.gff");
      System("mv chunk-preds.gff bulge$type-predictions-var$var.gff");
      addScores("bulge$type-predictions-var$var.gff",$sens,$snr);
      System("match_coordinates1.pl bulge$type-predictions-var$var.gff chunks-assembled bulge$type-predictions-var$var-genomic.gff");
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
`cut -f 2 $outfile-map.gff | sort | uniq -c > alignment_stats_$outfile-map`;

print "
# -------------------------------------------
# DONE.  OUTPUT IS IN: $outfile
# -------------------------------------------
";

print "
# -------------------------------------------
#  MOVING FILES TO DESIRED LOCATIONS CLEAN_UP
#  -------------------------------------------
";

#`mv *.gff $dir `;
#`mv chunks chunks-targetscan chunks-assembled $dir`;
#`mv out.maf miRNA_txt_cs miRNA_txt seeds.fasta PARCLIP.hmm $dir`;
#`mv file_defid.txt out_enst.txt targetscan_60_output.txt UTRs_median_BLs_bins.my_output.txt targetscan_60_output.sort.txt targetscan_60_output.BL_PCT.my_output.txt $dir`;
#`mv transcript-ids.txt targetscan_60_context_scores_output.txt tgscn_combined_out.txt $dir`;

`rm  transcript-ids.txt tgscn_combined_out.txt UTRs_median_BLs_bins.my_output.txt *.clusters.gff file_defid.txt out_enst.txt miRNA_txt* out.maf seeds.fasta *.tmp targetscan*`;
`rm -rf chunks* filename filename1 filename2 filename3 test.gff chunk-preds.gff *.hmm *.pl *.pm metamodel.txt submodels.txt cons.schema nocons.schema baum-welch classify extract-chain extract-motifs get-likelihood hmm-edit hmm-extract-state install-chain kmeans model-combiner parse random-HMM sample sub xgraph`;

System("date");


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

