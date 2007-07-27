# TAV and codon names for debugging:
# http://shadytrees.pastebin.ca/raw/537444 http://shadytrees.pastebin.ca/raw/537676
package Fabio;
use strict; use warnings;
use Smooth;
use List::Util qw(min max);

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
            my $codon = lc shift @a;
            $everything->{$codon} = shift @a;
        }
    };
    
    # All tRNA availabilities for every codon except
    # the stop codons.
    my @values = ();
    for my $key (keys %$everything) {
        push @values, $everything->{$key} unless $key =~ /uag|uga|uaa/;
    }
    
    our $smallest = min @values;
    our $biggest = max @values;
    $everything;
}

use POSIX;
our ($biggest, $smallest);
sub loops {
    my $num = shift;
    return ($biggest/$smallest - floor($num/$smallest));
    return ($biggest/$smallest) - floor($num/$smallest);
}

# Glorified print
sub plump {
    my ($info, $names) = @_;

    # Convert to lowercase because that's how
    # I stored the %everything keys, see above.
    my $line = 'Travel = struct(';
    Smooth::webopen $names, sub {
        s/\s//g; $_ = lc $_;
        $line .= "'$_', ";
        if ($_ =~ /uga|uag|uaa/) { $line .= '1000, ' }
        else { $line .= loops($info->{$_}) . ', ' }
    };
    $line =~ s/, $/)/;
    print $line;
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

    fabio.pl http://shadytrees.pastebin.ca/raw/537444 http://shadytrees.pastebin.ca/raw/537676 > temp.txt
    
=over 20

=item B<fabio.pl>

I<list of tRNA availabilites> I<list of codons>

=back

=head1 DESCRIPTION

fabio.pl's primary accomplishment is in parsing the
list of tRNA availabilities. It then sorts them in the
codon order specified by the list of codons, running them through
a Perl version of C<nloopcalcify> in order to generate a C<Travel.mat>
for GSPtools to use. Copy the output into a Matlab variable and use
C<save> to save that variable to a Travel.mat file.

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