#!/usr/bin/perl
# Extracts ONLY the character sequence from manually-created file

use strict;
my $seq = ''; my @seq = ( );

# Input file
my $infile = $ARGV[0];
# Read from file
open(INFILE, "$infile") or die "Cannot open file \"$infile\" for reading!\n\n";
    @seq = <INFILE>;
close(INFILE);

# Clean-up
$seq = join('', @seq);
# Remove whitespace and line numbers from sequence
$seq =~ s/[\s0-9]//g;

print $seq;
exit;