# 537444, 537676
# Usage:
#   * Obtain a list of frequencies.
#   * Convert Codons.mat to a text file with each codon separated by a newline.
#   * Pass the frequencies list and the codons list into fabio.pl.
#   * E.g. `fabio.pl frequencies.txt names.txt`

use strict;
use warnings;
package Fabio;

sub frequencies {
    my $freq = shift @_;
    open(my $handle, $freq) or die 'Cannot open frequencies';
    
    my @everything = ();
    while(<$handle>) {
        s/\(.*?\)//g;
        my @a = split /\s+/;
        grep !/^\s*$/, @a;
        push @everything, @a;
    }
    
    scalar {@everything};
}

sub plump {
    my ($info, $names) = @_;
    my %info = %$info;

    open my $handle, $names or die 'Cannot open names file';
    while (<$handle>) {
        s/\s//g; tr/a-z/A-Z/;
        print $info{$_}/1000, $/;
    }
    return;
}

sub main {
    my ($freq, $names) = @ARGV;
    plump(frequencies($freq), $names);
}

main() if __FILE__ eq $0; 1;