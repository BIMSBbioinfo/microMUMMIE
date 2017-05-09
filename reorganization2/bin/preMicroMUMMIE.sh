#!/usr/bin/env bash


###################################
##   MODIFY PATHS ACCORDINGLY    ##
###################################

#Put matchMIRS.pl and filterGroupsforMUMMIE.pl in the scripts directory and indicate path below
scripts=/data/ohler/Neel/bioTools/

#Download miRNA gff from ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 and indicate path below
mirnaGFF=/data/ohler/Neel/Genomes/human/hg19/annotation/v20/hsa.gff3

#Download miRNA fasta from ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz and indicate path below
mirnaFasta=/data/ohler/Neel/Genomes/human/hg19/annotation/v20/mature.fa



###################################
## DO NOT MODIFY ANYTHING BELOW  ##
###################################

cwd=$(pwd)
gff=$(basename "$mirnaGFF");
species=${gff%.*};
grep -v  miRNA_primary_transcript $mirnaGFF > $cwd/$species.p.gff3;
$cwd/$species.p.gff3;
for file in $1
do
filename=$(basename "$file")
library=${filename%.*}
if [ ! -d "$cwd/$library.$2.$3" ]; then
  mkdir "$cwd/$library.$2.$3";
else
 rm $cwd/$library.$2.$3/*;
fi

#calculate microRNA 5p/3p counts;
grep miRNA "$cwd/$file" | awk 'BEGIN { FS = "," } ; { print $1"\t"$3-1"\t"$4"\t"$14"\t"$8"\t"$2 }' | intersectBed -s -wao -a stdin -b $cwd/$species.p.gff3 | sed 's/;/\t/g' | sed 's/Name=//g' | awk '$NF > 17'  | awk '{a[$17]+=$5}END{for (i in a) print i,a[i]}' | sort -rnk2 | head -$2 > $cwd/$library.$2.$3/$library.$2.count.txt
wait;

#retrieve sequences in multiple formats for miRNAs found
ruby -e 'first_line = true; while line = STDIN.gets; line.chomp!; if line =~ /^>/; puts unless first_line; print line[1..-1]; print " "; else; print line; end; first_line = false; end; puts' < $mirnaFasta  | awk 'BEGIN{FS=" "} {print $1 "\t" $6}' > $cwd/$library.$2.$3/tmp
wait;
perl $scripts/matchMIRS.pl $cwd/$library.$2.$3/tmp $cwd/$library.$2.$3/$library.$2.count.txt $cwd/$library.$2.$3/$library.$2;
wait;
rm $cwd/$library.$2.$3/tmp;

#filter groups based on Conversion Specificity
perl $scripts/filterGroupsforMUMMIE.pl $cwd/$library.groups.csv $3 | sed 's/Aligned to/Aligned_to/' > $cwd/$library.$2.$3/$library.groups.csv
ln -s "$cwd/$library.distribution"  "$cwd/$library.$2.$3"
rm $cwd/$species.p.gff3;
done;
exit 0;
