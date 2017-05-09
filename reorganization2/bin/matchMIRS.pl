#!/usr/bin/env perl
use strict;
use warnings;

# $miRseqs is a csv file:
# GeneID,target sequence
# $rnaifile is a csv file:
# miRNA or siRNA name,sequence ***must be 5' to 3'
# if you want to look at duplex have separate entries for guide and passenger.
# for each miRNA/siRNA strand will output the 7mer seed match count and whether
# the target sequence has no match, 1 7mer match "7mer", >1 7mer matches "7mer_multi", or 7mer1A

(my $miRseqs, my $miRNAs, my $library ) = @ARGV;

#my $library = $prefix;
#Save the miRNA IDs of top x miRNAs in PAR-CLIP data
open(USRFL, "<",$miRNAs) or die "can't open $miRNAs $!\n";
my @seeds = ();
while(<USRFL>) {
	chomp $_;
	my @line = split(/\s/, $_);
	my $miRNA_ID = $line[0];
	push (@seeds,$miRNA_ID);
}
close USRFL;


#Match miRNAs found and write output files
open(USRFL, "<",$miRseqs) or die "can't open $miRseqs $!\n";
open (tabFILE, ">$library.miRNA.txt");
open (fastaFILE, ">$library.miRNA.fa");
while(<USRFL>) {
	chomp $_;
	my @line = split(/\t/, $_);
	my $gene = $line[0];
	my $seq = $line[1];
#	print "$gene,$seq\n";
	foreach my $seed ( @seeds) {	
		if ( $seed eq $gene) {
		print tabFILE "$gene\t$seq\n";
		print fastaFILE ">$gene\n$seq\n";
#		print "$gene\t$seq\n";
			}
	}
}
close (tabFILE);
close (fastaFILE);
close USRFL;

