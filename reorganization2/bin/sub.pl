#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <filename> \"old-text\" \"new-text\"\n" unless @ARGV==3;
my ($filename,$oldText,$newText)=@ARGV;

system("sed -i 's/$oldText/$newText/g' $filename");


