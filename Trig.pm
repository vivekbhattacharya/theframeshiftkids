use strict; use warnings;
package Trig;
use constant {
    PI => 3.14159,
    PERIOD => 4,
};

# This allows an algebraic trick for more elegant code.
our $canoot = 2 * PI / PERIOD;
sub fsin {
    my $x = $canoot * shift;
    my $phase = (PI / 2) - $canoot;
    
    if ($x > 4 - PERIOD and  $x < 4 + PERIOD) {
        return sin($x / 2 + $phase);
    }
    0;
}

sub bsin {
    my $x = $canoot * shift;
    my $phase = (PI / 2) + $canoot;
    
    if ($x > -4 - PERIOD and  $x < -4 + PERIOD) {
        return sin($x / 2 + $phase);
    }
    0;
}

sub xcos {
    my $x = $canoot * shift;
    if ($x > -1 * PERIOD and  $x < PERIOD) {
        return cos($x / 2);
    }
    0;
}

1;