# Frequency example: http://shadytrees.pastebin.ca/raw/537444
# Codon names example: http://shadytrees.pastebin.ca/raw/537676
# Usage:
#   * Obtain a list of frequencies.
#   * Convert Codons.mat to a text file with each codon separated by a newline.
#   * Pass the frequencies list and the codons list into fabio.pl.
#   * E.g. `fabio.pl frequencies.txt names.txt`
#   * Alternatively, post them on the Internet as text files (pastebins),
#     and pass the URLs.
# Output:
#   * TAV values in the order given by `names.txt`

use strict;
use warnings;
use Smooth;
package Fabio;

# Builds a map between codon and frequency.
sub frequencies {
    my $freq = shift;
    my @everything = ();
    
    # Remove parenthetical numbers and parse.
    Smooth::webopen $freq, sub {
        s/\(.*?\)//g;
        my @a = split /\s+/;
        @a = grep !/^\s*$/, @a;
        # Take advantage of Perl's array-as-hash hack.
        push @everything, @a;
    };
    
    # This is how you convert an array
    # to a hash reference?
    scalar {@everything};
}

# Glorified print
sub plump {
    my ($info, $names) = @_;

    # Convert to uppercase because frequencies
    # come in uppercase format.
    Smooth::webopen $names, sub {
        s/\s//g; tr/a-z/A-Z/;
        # Convert hash reference to hash first,
        # for those of you not fluent in sigilism.
        print ${%$info}{$_}/1000, $/;
    };
}

sub main {
    my ($freq, $names) = @ARGV;
    plump(frequencies($freq), $names);
}

main() if __FILE__ eq $0; 1;