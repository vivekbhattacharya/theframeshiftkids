# See Freier (1986)
package Kidnap::Freier;
use strict; use warnings;
use File::Basename;
use lib dirname(__FILE__);

use Kidnap::Util;
use constant BIG_NUM => 30000;

sub new {
    bless {};
}

sub init_penalty { 3.4 }

# Freier 1986
sub internal {
    my ($self, $t5, $t3, $b3, $b5) = @_;
    return BIG_NUM unless
      Kidnap::Util::valid_pair($t5, $b3)
          and Kidnap::Util::valid_pair($t3, $b5);

    my %scores =
        (
         AU => {
                # Watson/Crick matches
                AU => -0.9, UA => -0.9, GC => -1.7, CG => -2.1,
                # G/U mismatches
                GU => -0.5, UG => -0.7
               },
         UA => {
                AU => -1.1, UA => -0.9, GC => -1.8, CG => -2.3,
                GU => -0.7, UG => -0.5
               },
         CG => {
                AU => -1.8, UA => -1.7, GC => -2.9, CG => -2.0,
                GU => -1.5, UG => -1.5
               },
         GC => {
                AU => -2.3, UA => -2.1, GC => -2.9, CG => -3.4,
                GU => -1.3, UG => -1.9
               },
         GU => {
                AU => -0.5, UA => -0.7, GC => -1.5, CG => -1.9,
                GU => -0.5, UG => -0.5
               },
         UG => {
                AU => -0.7, UA => -0.5, GC => -1.5, CG => -1.3,
                GU => -0.6, UG => -0.5
               },
        );
    my $score = $scores{$t5 . $b3}->{$t3 . $b5};
    return $score if $score;
    die "Freier#internal(): Unable to pair top $t5$t3 with bottom $b3$b5";
}

sub terminal {
    # Terminal pairs must be at the very ends of both strands.
    my ($self, $t5, $t3, $b3, $b5, $left_side) = @_;

    # Flip bases if we are looking at the left side of the helix.
    if ($left_side) {
        ($t5, $b5) = ($b5, $t5);
        ($t3, $b3) = ($b3, $t3);
    }

    # We only need the first pair to match, provided it
    # is not a G/U mismatch.
    return BIG_NUM unless Kidnap::Util::wc_pair($t5, $b3);
    my %scores =
        (
         AU => {
                # Watson/Crick matches and G/U mismatches
                AU => -0.9, UA => -0.9, GC => -1.7, CG => -2.1,
                GU => -0.9, UG => -0.9,
                # Mismatches
                AA => -0.8, CC => -0.7, GG => -1.0, UU => -0.8,
                AC => -1.0, AG => -1.0, CA => -0.7, GA => -0.8,
                UC => -0.8, CU => -0.7,
               },
         UA => {
                AU => -1.1, UA => -0.9, GC => -1.8, CG => -2.3,
                GU => -0.9, UG => -1.0,
                AA => -1.0, CC => -0.6, GG => -1.2, UU => -0.5,
                AC => -0.8, AG => -1.1, CA => -0.7, GA => -1.1,
                UC => -0.6, CU => -0.5,
               },
         CG => {
                AU => -1.8, UA => -1.7, GC => -2.0, CG => -2.9,
                GU => -1.6, UG => -1.9,
                AA => -1.9, CC => -1.1, GG => -1.9, UU => -1.2,
                AC => -2.0, AG => -1.9, CA => -1.0, GA => -1.0,
                UC => -1.5, CU => -0.8,
               },
         GC => {
                AU => -2.3, UA => -2.1, GC => -2.9, CG => -3.4,
                GU => -1.4, UG => -2.3,
                AA => -1.1, CC => -0.6, GG => -1.4, UU => -0.7,
                AC => -1.7, AG => -1.3, CA => -1.1, GA => -1.6,
                UC => -0.8, CU => -0.5,
               },
        );
    my $score = $scores{$t5 . $b3}->{$t3 . $b5};
    return $score if $score;
    die "Freier#internal(): Unable to pair top $t5$t3 with bottom $b3$b5";
}

use Pod::Usage;
pod2usage(-verbose => 3, -output => \*STDOUT) if (__FILE__ eq $0);

1;

__END__

=head1 NAME

Kidnap::Freier

=head1 SYNOPSIS

    use Kidnap::Bind;
    my $o = Kidnap::Bind->new(Kidnap::Freier->new(37 + 273.15))
    $o->bind($seq, $rRNA);
    my $energy = $o->total_free_energy;

=head1 DESCRIPTION

Kidnap::Freier compiles values from Freier 1984 regarding the energy
of doublet bindings, both internal and doublet. Though not the most
current values, they are the only ones that work with GSPtools for now
due to all the work built upon them.

This module implements the API needed for a free energy module in
Kidnap.

=head1 METHODS

=over

=item $o->init_penalty()

Initial penalty in free energy, tacked on to the sum total in
Kidnap::Bind

=item $o->terminal_doublet()

=item $o->internal_doublet()

=back

=head1 COPYRIGHT AND LICENSE

See C<_Copying>.

=head1 SEE ALSO

Freier et al., PNAS (1986) Vol. 83, pp. 9373-9377

=cut
