library='test_lib'

preMicroMUMMIE.sh hsa.gff3 mature.fa AGO2_Z1.readcsv AGO2_Z1.groups.csv AGO2_Z1.distribution 100 1 $library
microMUMMIE.pl ./$library.100.1/$library.100.miRNA.txt Genome.2bit AGO2_Z1.groups.csv AGO2_Z1.distribution AGO2_Z1.clusters.csv  $library $library.merged.gff  0 ../reference/orderedtc_scUTRs.txt ./
