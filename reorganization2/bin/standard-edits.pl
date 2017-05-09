#!/usr/bin/perl
use strict;

system("hmm-edit PARCLIP.hmm DTRK phastcons VAR all 0 0.01");

# bin/hmm-edit PARCLIP.hmm TRK targetscan
#bin/hmm-edit PARCLIP.hmm MEAN all -- 1 0
#bin/hmm-edit PARCLIP.hmm VAR  all -- 1 3
#bin/hmm-edit PARCLIP.hmm MEAN 3   -- 1 3

