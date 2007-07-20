# TAV and codon names for debugging:
# http://shadytrees.pastebin.ca/raw/537444 http://shadytrees.pastebin.ca/raw/537676
package Fabio;
use strict;
use warnings;
use Smooth;

# Builds a map between codon and frequency.
sub frequencies {
    my $freq = shift;
    my $everything = {};
    
    # Remove parenthetical numbers and parse.
    use Data::Dumper;
    Smooth::webopen $freq, sub {
        s/\(.*?\)//g; s/^\s+//; s/\s+$//;
        my @a = split /\s+/;
        while (@a) {
            my $codon = shift @a;
            $everything->{$codon} = shift @a;
        }
    };
    $everything;
}

# Glorified print
sub plump {
    my ($info, $names) = @_;

    # Convert to uppercase because frequencies
    # come in uppercase format.
    Smooth::webopen $names, sub {
        s/\s//g; tr/a-z/A-Z/;
        print $info->{$_}/1000, $/;
    };
}

if (__FILE__ eq $0) {
    Smooth::helpcheck();
    my ($freq, $names) = @ARGV;
    plump(frequencies($freq), $names);
}
1;

__END__

=head1 NAME

fabio.pl (web-enabled)

=head1 SYNOPSIS

    fabio.pl http://shadytrees.pastebin.ca/raw/537444 http://shadytrees.pastebin.ca/raw/537676
    
=over 20

=item B<scan_brightly.pl>

I<list of tRNA availabilites> I<list of codons>

=back

=head1 DESCRIPTION

fabio.pl's primary accomplishment is in parsing the
list of tRNA availabilities. It then sorts them in the
codon order specified by the list of codons, making
this a suitable Matlab vector to run through C<nloopcalcify>
in order to generate a C<Travel.mat> for GSPtools to use.

=head1 OPTIONS AND ARGUMENTS

=over

=item I<list of availabilites>

You can obtain these from the Codon Usage Database[1]. Just
highlight and copy a portion of the text in the format seen
here[2].

=over

=item [1]

http://www.kazusa.or.jp/codon/

=item [2]

http://theframeshiftkids.googlecode.com/files/TAV.txt

=back

=item I<list of codons>

This must be a text file, each line containing a codon.
It must contain all the codons. This is the order in which
fabio.pl will output the TAV values.

=back

=cut