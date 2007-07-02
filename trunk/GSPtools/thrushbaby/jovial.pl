use strict; use warnings;
package Jovial;
use Smooth qw(prot2codon codon2prot);

sub shuffle {
    my ($seq, $times) = @_;
    my @codons = split /\s+/, $seq;
    my @proteins = codon2prot @codons;
    
    my @shuffles = ();
    for (my $i = 0; $i < $times; $i++) {
        push @shuffles, [prot2codon @proteins];
    }
    @shuffles;
}

if ($0 eq __FILE__) {
    # Schwartzian transform
    my @data = map { $_->[1] }
        sort { $a->[0] <=> $b->[0] }
        map {
            my ($codon, $loc, $times) = split '_';
            [$times, [$codon, $loc]];
        } @ARGV;
    use Data::Dumper;
    print Dumper(\@data);
}
1;