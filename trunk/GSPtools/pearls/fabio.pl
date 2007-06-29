# Frequency example: http://shadytrees.pastebin.ca/raw/537444
# Codon names example: http://shadytrees.pastebin.ca/raw/537676

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

our $help = <<END;
NAME
    fabio.pl
    (web-enabled)

USAGE
    perl fabio.pl [path to frequencies] [path to codons]
    See http://shadytrees.pastebin.ca/raw/537444 and
    http://shadytrees.pastebin.ca/raw/537676 for examples
    of the formats for these two files.
    
    This returns a list of frequencies that is suitable
    for creating TAV.mats. The list will follow the same
    order as the list of codons.
END
sub main {
    if (!@ARGV or $ARGV[0] eq '--help') { print $help; return }
    my ($freq, $names) = @ARGV;
    plump(frequencies($freq), $names);
}

main() if __FILE__ eq $0; 1;