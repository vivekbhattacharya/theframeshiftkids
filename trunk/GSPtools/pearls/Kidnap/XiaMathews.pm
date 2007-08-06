# See Xia (1998) (for Watson-Crick pairs) and Mathews (1999) (for G/U mismatches)
use File::Basename;
use lib dirname(__FILE__);
use strict; use warnings; use Util;

package Kidnap::XiaMathews;
use constant BIG_NUM => 30000;

sub new {
	my ($class, $temp) = @_;
    my $self = {InitPenalty => 3.4, Temp => $temp};
    bless($self, $class);
    
    my ($dH, $dS) = (3.61, -1.5);
    $self->{InitPenalty} = $dH - $temp*($dS/1000);
    return $self;
}

sub init_penalty { shift->{InitPenalty} }

sub internal {
	my $self = shift;
    my ($t5, $t3, $b3, $b5, $t55, $t33, $b33, $b55) = @_;
    unless (Util::valid_pair($t5, $b3) && Util::valid_pair($t3, $b5)) {
		return BIG_NUM;
    }
    
	my %scores = (
		AU => {
			AU => [-6.82, -19.0], UA => [-9.38, -26.7],
			GC => [-10.48, -27.1], CG => [-11.40, -29.5],
			GU => [-3.21, -8.6], UG => [-8.81, -24.0],
		}, UA => {
			AU => [-7.69, -20.5], UA => [-6.82, -19.0],
			GC => [-10.44, -26.9], CG => [-12.44, -32.5],
			GU => [-6.99, -19.3], UG => [-12.83, -37.3],
		}, CG => {
			AU => [-10.44, -26.9], UA => [-10.48, -27.1],
			GC => [-10.64, -26.7], CG => [-13.39, -32.7],
			GU => [-5.61, -13.5], UG => [-12.11, -32.2],
		}, GC => {
			AU => [-12.44, -32.5], UA => [-11.40, -29.5],
			GC => [-13.39, -32.7], CG => [-14.88, -36.9],
			GU => [-8.33, -21.9], UG => [-12.59, -32.5],
		}, GU => {
			GU => [-13.47, -41.82], UG => [-14.59, -51.2],
			UA => [-8.81, -24.0], AU => [-12.83, -37.3],
			GC => [-12.11, -32.3], CG => [-12.59, -32.5],
		}, UG => {
			GU => [-9.26, -30.8], UG => [-13.47, -41.82],
			UA => [-3.21, -8.6], AU => [-6.99, -19.3],
			GC => [-5.61, -13.5], CG => [-8.33, -21.9],
		},
	);

	# Mathews (1999) scores the entire GGUC/CUGG unit with
	# dH = -30.80 and dS = -86.0. The values below
	# are via
	# 	dH = -30.80-(2*-8.33) where -8.33 is dH for GG/CU
	# 	dS = -86-(2*-21.9) where -21.9 is dS for GG/CU
	my ($dH, $dS) = @{$scores{$t5 . $b3}->{$t3 . $b5}};
	if (defined $b55 and join('', @_) eq 'GCCG') {
		($dH, $dS) = (-14.14, -42.2);
	}
    return $dH - $self->{Temp}*($dS/1000);
}

sub terminal {
	my $self = shift;
    my $doublet_score = $self->internal(@_);

	my ($t5, $t3, $b3, $b5, $left_side) = @_;
	my @chesthair = $left_side ? ($t5, $b3) : ($t3, $b5);
	my $terminal_penalty = $self->terminal_pair(@chesthair);

    return $doublet_score + $terminal_penalty;
}

## This scores a pair of RNA residues based on theoretical binding
## ability, assuming the pair exists at the terminal end of a helix
## structure. Terminal pairs must be at the ends of both strands.
sub terminal_pair {
    my ($self, $t, $b) = @_;

    my $x = $t . $b;
	if (grep /$x/, qw/AU UA GU UG/) {
		my ($dH, $dS) = (3.72, 10.5);
        return $dH - $self->{Temp}*($dS/1000);
    }
    else { return 0; }
}
1;