#!/usr/bin/perl
#Program written to print genomic coordinates from UTR file in the output gff file September 27 2013 Samta Malhotra
use ProgramName;
use strict;

my $name=ProgramName::get();
die "$name <in.gff> <directory fastb files> <out.gff>\n" unless @ARGV==3;
my ($infile,$dir,$outfile)=@ARGV;
my @data;
my @ensemblid;
my @files=`ls $dir`;
my $n=@files; # length of all fastb files
open(IN,$infile) || die "can't open file: $infile\n";
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
while(<IN>) {
  chomp;
  my @fields=split(/\s+/,$_);
  my $substrate=$fields[0];
  my @value = split(/_/,$substrate);
  my $GeneID = $value[0];
  my $begin=$fields[3];
  my $end=$fields[4];
  my $score=$fields[5];
  my $strand=$fields[6];
  my $dot=$fields[7];
  my $extra=$fields[8];
  #print "The gene is=$GeneID,\n";
  for(my $j=0 ; $j<$n ; ++$j) {  #as long as $j is smaller $n (number of fastb files)

    # define variable read and defined from fastb file
    my @exonstart2;
    my @exonend2;
    my @split;
    my $all;
    my $size;
    my $diff;
    my $numexons;
    my $si,
    my $newbegin;
    my $newend;
    my $le;
    my $mystartvalue;
    my $myendvalue;
    my $exonstart2;
    my $exonend2;
    my $file=$files[$j]; # take each fastb file one by one; by $j
    chomp ($file);
    #print "$file\t$substrate.fastb\n";
    if ("$file" eq "$substrate.fastb") {
    #if ($file =~ m/(^$substrate).fastb/) { # if fastb file handle matches the gene name
      #print "$j\n";
      my $name1=$substrate;
      my $sscore = $score;
      open(FILE,"$dir/$substrate.fastb") || die "can't open $dir/$substrate.fastb\n"; # open fastb file
      while(<FILE>) {
        #if($_ =~ m/transcriptID=(ENST\d+)/) ; hard coded, only works for ENSEMBL transcript IDs
        if($_ =~ m/transcriptID=(FBtr\d+)/) { # find line with TranscriptID
          my $transcriptid = $1;  # define transcript ID
          if ($_ =~ m/numexons=(\d+)/){
            $numexons = $1; # define exon numer
            if ($_ =~ m/chr=(.*)/){
              my $chrm = $1; # define chromosome
              my @splitch = split(/\s+/,$chrm);
              my $chr = $splitch[0];
              if ($_ =~ m/strand=(.*)/){
                $all = $1;
                @split = split(/\s+/,$all);
                my $sstrand = $split[0];
                my $sum=0;
                $size = scalar(@split);
                my $size1 = $size-1;
                my @exonstart1;
                my @exonend1;
                for (my $s=1;$s<$size1;$s=$s+2){
                  push (@exonstart1,$split[$s]);
                  $s=$s+1;
                  push (@exonend1,$split[$s]);
                  $s=$s-1;
                }
                my @exonStart;
                my @exonEnd;
                my @exonstart;
                my @exonend;
                my $sum=0;
                my $newdiff=0;
                $le = scalar(@exonend1);
                for (my $a=0;$a<$le;$a++){
                  @exonstart2 = split(/=/,$exonstart1[$a]);
                  @exonend2 = split(/=/,$exonend1[$a]);
                  push (@exonStart,$exonstart2[1]);
                  push (@exonEnd,$exonend2[1]);
                  @exonstart = sort @exonStart;
                  @exonend = sort @exonEnd;
                }
                for (my $b=0;$b<$le;$b++){
                  my $se=$le-1;
                  if ($sstrand eq "+"){
                    $diff = $exonend[$b]-$exonstart[$b];
                    $newdiff = $newdiff+$diff;
                  }
                  if ($sstrand eq "-"){
                    my $sel=$se-$b;
                    $diff = $exonend[$sel]-$exonstart[$sel];
                    $newdiff = $newdiff+$diff;
                  }
                  if ($begin>$newdiff){
                    $sum = $sum+$diff;
                    my $newBegin = $begin-$sum;
                    my  $newEnd = $end-$sum;
                    $newbegin = $newBegin;
                    $newend = $newEnd;
                  }
                  else{
                    my ($chunkBegin,$chunkEnd)=($1,$2);
                    if ($le<=1){
                      my $lenn=$le-1;
                      if($sstrand eq "+") {
                        $chunkBegin=$exonstart[0]+$begin;
                        $chunkEnd=$exonstart[0]+$end;
                      }
                      if ($sstrand eq "-") {
                        # $mystartvalue = $exonend[$lenn]-$begin;
                        # $myendvalue=$exonend[$lenn]-$end;
                        $mystartvalue = $exonend[$lenn]-$begin+1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                        $myendvalue=$exonend[$lenn]-$end+1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                        $chunkBegin=$myendvalue;
                        $chunkEnd=$mystartvalue;
                      }

                      #print OUT "chr$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$lenn];prediction_start=$begin;prediction_end=$end;$extra\n"; last;

                      print OUT "$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$lenn];prediction_start=$begin;prediction_end=$end;$extra\n"; last;
                    }

                    if ($le>1){
                      if (($b==0) and ($begin < $diff)){
                        my $len=$le-1;
                        if($sstrand eq "+") {

                          $chunkBegin=$exonstart[0]+$begin;
                          $chunkEnd=$exonstart[0]+$end;

                          if ($chunkEnd > $exonend[0]){
                            my $s=$b+1;
                            my $difp=$chunkEnd-$exonend[0];
                            my $valtoo = $exonstart[$s];
                            my $valll = $valtoo+$difp;
                            $chunkEnd=$exonend[0];
                            #print OUT "chr$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$len];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";

                            print OUT "$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$len];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";
                          }
                        }
                        if ($sstrand eq "-") {

                          # $mystartvalue = $exonend[$len]-$begin;
                          # $myendvalue=$exonend[$len]-$end;
                          $mystartvalue = $exonend[$len]-$begin+1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          $myendvalue=$exonend[$len]-$end+1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          if ($myendvalue < $exonstart[$len]){
                            my $v=$len-1;
                            my $dif=$exonstart[$len]-$myendvalue;
                            my $valto = $exonend[$v];
                            my $vall = $valto-$dif;
                            $myendvalue=$exonstart[$len];
                            #print OUT "chr$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction$extra\n";

                            print OUT "$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction$extra\n";
                          }
                          $chunkBegin=$myendvalue;
                          $chunkEnd=$mystartvalue;
                        }

                        #print OUT "chr$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$len];prediction_start=$begin;prediction_end=$end;$extra\n"; last;

                        print OUT "$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$len];prediction_start=$begin;prediction_end=$end;$extra\n"; last;
                      }

                      elsif (($b>0) and ($begin>$newdiff)){
                        my $se=$le-1;
                        if($sstrand eq "+") {
                          my $a=$b+1;
                          $chunkBegin=$exonstart[$b]+$newbegin;
                          $chunkEnd=$exonstart[$b]+$newend;
                          # print "my end is=$chunkEnd,\n";
                          if ($chunkEnd>$exonend[$b]){
                              my $difp=$chunkEnd-$exonend[$b];
                              my $valtoo = $exonstart[$a];
                              my $valll = $valtoo+$difp;
                              $chunkEnd=$exonend[$b];
                              #print OUT "chr$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";

                            print OUT "$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";
                          }
                        }
                        if ($sstrand eq "-") {
                          my $sel=$se-$b;
                          # $mystartvalue = $exonend[$sel]-$newbegin;
                          # $myendvalue=$exonend[$sel]-$newend;
                          $mystartvalue = $exonend[$sel]-$newbegin+1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          $myendvalue=$exonend[$sel]-$newend+1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          if ($myendvalue<$exonstart[$sel]){
                            my $v=$sel-1;
                            my $dif=$exonstart[$sel]-$myendvalue;
                            my $valto = $exonend[$v];
                            my $vall = $valto-$dif;
                            #print OUT "chr$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";

                            print OUT "$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";
                            $myendvalue=$exonstart[$sel];
                          }
                          $chunkBegin=$myendvalue;
                          $chunkEnd=$mystartvalue;
                        }
                        #print OUT "chr$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;$extra\n"; last;

                        print OUT "$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;$extra\n"; last;
                      }

                      elsif (($b>0) and ($begin<$newdiff)){
                        my $se=$le-1;
                        if($sstrand eq "+") {
                            my $a=$b+1;
                          $chunkBegin=$exonstart[$b]+$newbegin;
                          $chunkEnd=$exonstart[$b]+$newend;
                          if ($chunkEnd>$exonend[$b]){
                            my $difp=$chunkEnd-$exonend[$b];
                            my $valtoo = $exonstart[$a];
                            my $valll = $valtoo+$difp;
                            $chunkEnd=$exonend[$b];
                            #print OUT "chr$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";

                            print OUT "$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";
                          }
                        }
                        if ($sstrand eq "-") {
                          my $sel=$se-$b;
                          #$mystartvalue = $exonend[$sel]-$newbegin;
                          #$myendvalue=$exonend[$sel]-$newend;
                          $mystartvalue = $exonend[$sel]-$newbegin + 1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          $myendvalue=$exonend[$sel]-$newend + 1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          if ($myendvalue<$exonstart[$sel]){
                            my $v=$sel-1;
                            my $dif=$exonstart[$sel]-$myendvalue;
                            my $valto = $exonend[$v];
                            my $vall = $valto-$dif;
                            #print OUT "chr$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";

                            print OUT "$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";
                            $myendvalue=$exonstart[$sel];
                          }
                          $chunkBegin=$myendvalue;
                          $chunkEnd=$mystartvalue;
                        }
                        #print OUT "chr$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;$extra\n"; last;

                        print OUT "$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;$extra\n"; last;
                      }


                      elsif (($b>0) and ($begin==$newdiff)){
                        my $se=$le-1;
                        if($sstrand eq "+") {
                          my $a=$b+1;
                          $chunkBegin=$exonstart[$b]+$newbegin;
                          $chunkEnd=$exonstart[$b]+$newend;
                          if ($chunkEnd>$exonend[$b]){
                            my $difp=$chunkEnd-$exonend[$b]-1;
                            my $valtoo = $exonstart[$a];
                            my $valll = $valtoo+$difp;
                            $chunkEnd=$exonend[$b];
                            #print OUT "chr$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";

                            print OUT "$chr\tBinding\tsite\t$valtoo\t$valll\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";
                          }
                        }
                        if ($sstrand eq "-") {
                          my $sel=$se-$b;
                          # $mystartvalue = $exonend[$sel]-$newbegin;
                          # $myendvalue=$exonend[$sel]-$newend;
                          $mystartvalue = $exonend[$sel]-$newbegin + 1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          $myendvalue=$exonend[$sel]-$newend + 1; #This was added by harm; becasue predictions on the "-" strand were shifted by -1
                          if ($myendvalue<$exonstart[$sel]){
                            my $v=$sel-1;
                            my $dif=$exonstart[$sel]-$myendvalue;
                            my $valto = $exonend[$v];
                            my $vall = $valto-$dif;
                            #print OUT "chr$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";

                            print OUT "$chr\tBinding\tsite\t$vall\t$valto\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;spliced_prediction;$extra\n";
                            $myendvalue=$exonstart[$sel];
                          }
                          $chunkBegin=$myendvalue;
                          $chunkEnd=$mystartvalue;
                        }
                        #print OUT "chr$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;$extra\n"; last;

                        print OUT "$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$se];prediction_start=$begin;prediction_end=$end;$extra\n"; last;
                      }
                    }
                    #print "Exonstart = $exonstart[$b], Predictionstart = $chunkBegin, Predictionend = $chunkEnd, Fileprediction = $begin, Fileend = $end, newpredictionstart =$newbegin, newpredictionend = $newend\n";
                    #print OUT "chr$chr\tBinding\tsite\t$chunkBegin\t$chunkEnd\t$score\t$sstrand\t$dot\tgene=$substrate;transcriptid=$transcriptid;utr_start=$exonstart[0];utr_end=$exonend[$b];prediction_start=$begin;prediction_end=$end;$extra\n";
                  }
                }
              }
            }
          }
        }
      }
      close(FILE);
    }
    else {$file=$files[$j+1];}
  }
}
close (IN);
close (OUT);
close (UTR);
