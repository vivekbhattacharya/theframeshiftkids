package Util;
use strict; use warnings;

## Returns true if two characters form a valid helix pair,
## i.e. a standard Watson-Crick pair or a G-U mismatch.
sub valid_pair {
    my ($a, $b) = @_;
	my $ab = $a . $b;
	scalar grep /$ab/, qw/GC CG AU UA GU UG/;
}


## Returns true if two characters form a valid W-C pair
sub wc_pair {
    my ($a, $b) = @_;
	my $ab = $a . $b;
	scalar grep /$ab/, qw/GC CG AU UA/;
}

sub dna2rna {
    my $seq = shift;
    $seq =~ tr/tT/uU/;
	return $seq;
}
1;