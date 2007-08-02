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
    
    my ($overaged, $type, $ant) = (0, 0, '');
    my $wait; my $age_limit = 150 + 100 * rand;
    for ($wait = 1; $wait < $age_limit + 1; $wait++) {
        my ($a, $b, $c) = map { $_->update($x0 - 2 * $shift) } ($back, $here, $next);
        
        # Careful with the single letters
        my $reloop = 1 * (1 - ($a + $b + $c));
        if ($reloop < $a or $reloop < $b or $reloop < $c) {
            my $r = rand;
            if ($r < $b) {
                $overaged = $here->is_stop;
                last;
            } elsif ($r < $b + $c) {
                $overaged = $next->is_stop;
                ($type, $ant) = (1, $codon);
                last;
            } elsif ($r < $a + $b + $c) {
                $overaged = $back->is_stop;
                ($type, $ant) = (2, $codon);
                last;
            }
        }
        
        $x0 += $secret->($x0);
    }
    
    $overaged = 2 if $wait > $age_limit;
    print "{ $x0 $wait $type '$ant' $overaged }", $/;
}
1;