use strict; use warnings;
package Codon;
use constant {
    POWER => 10,
    STOP => 1000,
};

sub new {
    my $class = shift;
    
    my ($loops, $weigh) = @_;
    my $self = {
        loops => $loops, weigh => $weigh,
        fail => 1,
    };
    bless $self, $class;
    return $self;
}

sub update {
    my ($self, $window) = @_;

    $window = $self->{weigh}->($window) ** POWER;
    $self->{fail} *= 1 -  $window / $self->{loops};
    return 1 - $self->{fail};
}

sub is_stop {
    my $self = shift;
    if ($self->{loops} == STOP) { return 1 }
    return 0;
}
1;