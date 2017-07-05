#!/usr/bin/env bash


###################################
##   MODIFY PATHS ACCORDINGLY    ##
###################################

#Download miRNA gff from ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 and indicate path below
mirnaGFF=$1

#Download miRNA fasta from ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz and indicate path below
mirnaFasta=$2

readcsv=$3

groups_csv=$4

distribution=$5

top=$6

specif=$7

library=$8

output=$library.$top.$specif

###################################
## DO NOT MODIFY ANYTHING BELOW  ##
###################################

cwd=$(pwd)

if echo "$mirnaGFF" | grep -q "gff3"; then
    grep  miRNA $mirnaGFF | grep -v miRNA_primary_transcript > $cwd/tmp.p.gff3;
fi

if [ ! -d "$cwd/$output" ]; then
  mkdir "$cwd/$output";
else
 rm $cwd/$output/*;
fi

#calculate microRNA 5p/3p counts;
grep miRNA "$cwd/$readcsv" | awk 'BEGIN { FS = "," } ; { print $1"\t"$3-1"\t"$4"\t"$14"\t"$8"\t"$2 }' | intersectBed -s -wao -a stdin -b $cwd/tmp.p.gff3 | sed 's/;/\t/g' | sed 's/Name=//g' | awk '$NF > 17'  | awk '{a[$17]+=$5}END{for (i in a) print i,a[i]}' | sort -rnk2 | head -$top > $cwd/$output/$library.$top.count.txt
wait;

#retrieve sequences in multiple formats for miRNAs found
python -c """
seqs = open('$mirnaFasta').readlines()

for i in range(0, len(seqs), 2):
    seq_id = seqs[i].split()[0][1:]
    seq = seqs[i+1].strip()
    print('%s\t%s' % (seq_id, seq))
""" > $cwd/$output/tmp
wait;
matchMIRS.pl $cwd/$output/tmp $cwd/$output/$library.$top.count.txt $cwd/$output/$library.$top;
wait;
rm $cwd/$output/tmp;

#filter groups based on Conversion Specificity
filterGroupsforMUMMIE.pl $cwd/$groups_csv $specif | sed 's/Aligned to/Aligned_to/' > $cwd/$output/$library.groups.csv
ln -s "$cwd/$distribution"  "$cwd/$output/$library.distribution"
rm $cwd/tmp.p.gff3;

exit 0;
