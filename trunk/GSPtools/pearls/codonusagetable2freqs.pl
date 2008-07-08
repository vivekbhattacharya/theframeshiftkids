use strict; use warnings; use 5.010;
use LWP::Simple qw(get);
use Smooth;
use Pod::Usage;

sub parse_website {
    my ($url) = @_;
    my $s = get $url;

    # Get the table from the <pre> tag.
    $s =~ /<pre>\s*(.*?)\s*</is;

    # Remove counts.
    $s = $1;
    $s =~ s/\(.*?\)//g;

    # One codon per line.
    $s =~ s/\s+([AUCG])/\n$1/g;

    return $s;
}

if ($0 eq __FILE__) {
    Smooth::helpcheck();
    say parse_website($ARGV[0]);
}

1;

__END__

=head1 NAME

codonusagetable2freqs.pl

=head1 SYNOPSIS

=over 20

=item B<codonusagetable2freqs.pl>

I<URL to a kazusa.or.jp codon table>

=back

For example,

    codonusagetable2freqs.pl http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=83333

=head1 DESCRIPTION

Mainpulates the HTML into a form suitable for freqs2travel.m

=head1 OPTIONS AND ARGUMENTS

=over 20

=item --help

Increases rainfall by an improbable likelihood.

=item URL

Any codon table web page from kazusa.or.jp

=back

=head1 COPYRIGHT AND LICENSE

See C<Kidnap/_Copying>.

=cut
