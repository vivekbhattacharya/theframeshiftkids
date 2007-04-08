#!/usr/bin/perl
# Extracts ONLY the character sequence from manually-created file

use strict;
my $file = shift @ARGV;
open(F, $file) or die "Cannot open file \"$file\" for reading!\n\n";
        $_ = join '', <F> and s/[\s0-9]//g and print;
close(F);