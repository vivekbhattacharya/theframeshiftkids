#!/usr/bin/env perl

use warnings; use strict; use 5.008;
use File::Basename;
use lib dirname(__FILE__);

use Getopt::Std;
use Smuggle::Bind;
use Smooth;
use Cache;

sub dna2rna {
    my $seq = shift;
    $seq =~ tr/tT/uU/;
    return $seq;
}

if ($0 eq __FILE__) {
    Smooth::helpcheck();

    our ($opt_t, $opt_n) = (37, 0);
    getopts('nt:');
    # Convert Celsius to Kelvin.
    $opt_t += 273.15;

    local @_ = @ARGV;
    my $rna  = shift
      or die 'scan_bindigo.pl: No RNA binding sequence given';
    my $seq  = shift
      or die 'scan_bindigo.pl: No sequence file given';

    # BINDIGO wants 5' to 3'.
    $rna = scalar reverse $rna;

    # 'type:bindigo' will separate the cached files from their
    # scan_brightly counterparts.
    $seq = Smooth::getseq $seq;

    my $cache;
    unless ($opt_n) {
        $cache = Cache->init('scan_brightly_cache');
        my $x  = $cache->try_get($rna, $seq, $opt_t, 'type:bindigo');
        print $x and exit if $x;
    }

    # If we're here, it's not cached.
    $rna = lc dna2rna($rna);
    $seq = lc dna2rna($seq);

    my $x = Smuggle::Bind::bind_str($rna, $seq);
    print($x);

    # Since it's not cached, cache it.
    $cache->new_store($x) unless $opt_n;
}
1;

__END__

=head1 NAME

scan_bindigo.pl

=head1 SYNOPSIS

    scan_bindigo.pl -t 273.15 auuccuccacuag prfB.fasta

=over 20

=item B<scan_brightly.pl>

[-n] [B<-p> I<parameter>] [B<-t> I<temperature (K)>] I<RNA sequence>
I<Gene file>

=back

=head1 DESCRIPTION

scan_brightly takes an RNA binding sequence (like the 3 prime tail on the 16s
rRNA) and calculates the free energy signals that determine how well it
will bind to the sequence (DNA or RNA) at every point, outputting it.
Each signal is separated by a newline, allowing Matlab to store the result
from a C<system> call as a vector.

=head1 OPTIONS AND ARGUMENTS

=over

=item --help

Increases rainfall by an improbable likelihood.

=item -p MODULE

This determines which set of parameters are used to simulate binding
between the nucleotides strands.  The default module is Freier (1986).
Alternatively, XiaMathews (1998) is more modern but does not work
GSPtools. scan_brightly will search C<@INC> for C<MODULE.pm>. Use Perl's
double-colon syntax, such as C<Kidnap::XiaMathews>.

=item -t NUM

This specifies the temperature, in Celsius, at which the binding occurs.
The default value is 37, human body temperature, for whatever reason.

=item -n

This turns off signal caching. It's useful for long runs (like the
entirety of the E. coli genome) when you fear for your hard drive.

=item RNA sequence

See below.

=item Gene file

This can be in a FASTA format or just the nucleotides.

=back

=head1 COMMON RNA SEQUENCES

Usually, the RNA sequences are 3' 16/18S tails. Here are several
written in the 3'->5' orientation.

=over 40

=item I<E. coli>

auuccuccacuag (13 bases)

=item I<Z. mays>

guuacuag (8 bases)

=item I<Arabidopsis thaliana>

guuacuag (8 bases)

=item I<Xenopus laevis>

auuacuag (8 bases)

=item I<Mus musculus>

auuacuag (8 bases)

=item I<D. melanogaster>

auuacuag (8 bases)

=item I<Homo sapiens>

auuacuag (8 bases)

=back

=head1 CAVEATS

The sequence data contained in FASTA file must be in a 5'->3'
orientation. The RNA sequence must be in a 3'->5' orientation.

=head1 AUTHOR

Originally by Joshua Starmer, scan_brightly has since then
been heavily modified by Hao Lian.

=head1 COPYRIGHT AND LICENSE

See C<Kidnap/_Copying>.

=cut
