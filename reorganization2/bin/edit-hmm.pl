#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <*.hmm>\n" unless @ARGV==1;
my ($outfile)=@ARGV;

system("MUMMIE/hmm-edit $outfile MEAN 0 0 1e-7   MEAN 0 1 0.531313");
system("MUMMIE/hmm-edit $outfile MEAN 1 0 0.5    MEAN 1 1 0.531313");
system("MUMMIE/hmm-edit $outfile MEAN 2 0 1.0    MEAN 2 1 0.531313");
system("MUMMIE/hmm-edit $outfile MEAN 3 0 0.5    MEAN 3 1 0.993844");
system("MUMMIE/hmm-edit $outfile COV ALL ALL ALL 0");
system("MUMMIE/hmm-edit $outfile VAR  0 0 0.01   VAR  0 1 1.10313");
system("MUMMIE/hmm-edit $outfile VAR  1 0 0.01   VAR  1 1 1.10313");
system("MUMMIE/hmm-edit $outfile VAR  2 0 0.01   VAR  2 1 1.10313");
system("MUMMIE/hmm-edit $outfile VAR  3 0 0.01   VAR  3 1 1.15067");

