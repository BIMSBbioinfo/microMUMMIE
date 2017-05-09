#!/usr/bin/env perl
use strict;
use warnings;

# Filter group file based on Conversion Specificity
(my $groups, my $CS_filter ) = @ARGV;

if ($CS_filter eq "none") {
	open(USRFL, "<",$groups) or die "can't open $groups $!\n";

	while(<USRFL>) {
	my @line = split(/,/, $_);
	print "$_";
}
close USRFL;
exit;
}


open(USRFL, "<",$groups) or die "can't open $groups $!\n";

	while(<USRFL>) {
	my @line = split(/,/, $_);

	my $CS = $line[13];
			
		if ( $_ =~ /ConversionSpecificity/) {
			print "$_";
				}
		elsif ($CS ne "NA" and $CS > $CS_filter) {
			print "$_";
		}
	
}

close USRFL;
