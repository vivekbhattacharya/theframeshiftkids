use strict; use warnings;
package Rebase;
use File::Basename;
use lib dirname(__FILE__);
use Freebase::Trig;
use Freebase::Codon;

if (__FILE__ eq $0) {
    my @loops = @ARGV[0..2];
    my ($codon, $x0, $shift) = @ARGV[3..5];
    my ($diff1, $diff2) = @ARGV[6..7];
    
    my $back = Codon->new(shift @loops, \&Trig::bsin);
    my $here = Codon->new(shift @loops, \&Trig::xcos);
    my $next = Codon->new(shift @loops, \&Trig::fsin);
    
    # This follows from phi_signal(1,k) = Dvec(k,2)
    # "A model for +1 frameshifts in eubacteria" by Ponnala, et al.
    my $species = -30 * (Trig::PI / 180);
    my $secret = sub {
        my $thigh = shift() * (Trig::PI / 3) - $species + $diff2;
        return -0.005 * $diff1 * sin($thigh);
    };
    
    my ($overaged, $ant, $termite) = (0, '', '');
    my $wait;
    for ($wait = 1; $wait < 150 + 100 * rand; $wait++) {
        my ($a, $b, $c) = map { $_->update($x0 + 2 * $shift) } ($back, $here, $next);
        
        my $reloop = 1 * (1 - ($a + $b + $c));
        if ($reloop < $a or $reloop < $b or $reloop < $c) {
            my $r = rand;
            if ($r < $a) {
                $overaged = $here->is_stop;
                last;
            } elsif ($r < $a + $b) {
                $overaged = $next->is_stop;
                $ant = $codon;
                last;
            } elsif ($r < $a + $b + $c) {
                $overaged = $back->is_stop;
                $termite = $codon;
                last;
            }
        }
        
        $x0 += $secret->($x0);
    }
    
    print $x0, $/;
    print $wait;
}
1;