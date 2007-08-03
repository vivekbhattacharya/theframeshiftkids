package Kidnap::Bind;
use warnings; use strict;
use File::Basename;
use lib dirname(__FILE__);

sub new {
	my $self = {};
    my ($class, $doublet) = @_;
	
    # Parameters module (pluggable)
    $self->{Doublet} = $doublet;
	
	bless($self, $class);
}

use Strand;
sub bind {
    my ($self, $seq, $rna) = @_;    
    my $doublet = $self->{Doublet};

	# I need at least two bases in each sequence in order to form
	# a structure between them.
    my $length = length $seq;
    if ($length != length $rna) {
		die 'Bind#bind: Unequal sequence lengths';
    } if ($length < 2) {
		die 'Bind#bind: Sequences are too short';
    }

    my $best_score = my $score = 0;
    my $helix_end = 0;
    
    my $x = Strand->new($seq, $length);
    my $y = Strand->new($rna, $length);
    for my $i (0 .. $length-2) {
        $x->update($i); $y->update($i);
		
		# Start each helix structure with a terminal doublet.
		if ($score == 0) {
			$score = $doublet->terminal($x->all, $y->all, 1);
		} else {
			# Some internal "doublets" need more context to
			# be scored correctly.  See the scoring for GU
			# in XiaMathews#internal for more details.	
			$score += $doublet->internal($x->all, $y->all, $x->context, $y->context);
		}
		
		if ($score > 0) { $score = 0 }
        elsif ($score < $best_score) {
			$best_score = $score;
			$helix_end = $i;
		}
    }

	# Swap out the last internal for a terminal doublet score
	# so that the helix starts and ends with a terminal doublet.
	# We have to do this now because we cannot know in advance
	# where this terminal doublet is going to be.
    if ($best_score < 0) {
		$x->update($helix_end); $y->update($helix_end);
		
		my $internal_score = $doublet->internal($x->all, $y->all);
		my $terminal_score = $doublet->terminal($x->all, $y->all, !1);
		$best_score -= $internal_score - $terminal_score;
	
        # We should check for "self symmetry." With the force_bind,
        # we are assuming that the strands are stretched out, and thus
        # not stuck to themselves forming their own hairpin loop. The 
        # "self symmetry" penalty arises from the individual strands
        # being stuck together in their own hairpins.
    }

    return $best_score + $doublet->init_penalty;
}

use Pod::Usage;
pod2usage(-verbose => 3, -output => \*STDOUT) if (__FILE__ eq $0);

1;

__END__

=head1 NAME

Kidnap::Bind 

=head1 SYNOPSIS

    use Kidnap::Bind;
    my $o = Kidnap::Bind->new(Kidnap::Freier->new(37 + 273.15))
    $o->bind($seq, $rRNA);
    my $energy = $o->total_free_energy;

=head1 DESCRIPTION

Kidnap::Bind calculates free energy by forcing an alignment
between rRNA and mRNA and scoring RNA residues.

=head1 TERMINOLOGY

A doublet consists of two pairs of RNA residues.

=head1 METHODS

=over

=item $o = Kidnap::Bind->new($parameters)

C<$parameters> is a parameter module such as C<Kidnap::Freier>
or C<Kidnap::XiaMathews> that implements C<internal_doublet>
and C<terminal_doublet> methods. These methods should compute
the score for pairs of RNA residues using numbers pulled from
any given paper such as Freier 1986.

=item $o->internal_doublet($t5, $t3, $b3, $b5)

C<$t5> and C<$t3> represent the top strand's 5' and 3'. C<$b3>
and C<$b5> represent the bottom strand's 3' and 5'.

This returns free energy score of a  doublet of RNA residue based on 
how well they would bind, assuming that the doublet is found inside a
helix structure. Thus, this only considers Watson/Crick paris plus G/U
pairs. 

This only scores matches or internal (to a helix structure, i.e.
not in a loop or bulge) G/U mismatches.  If you're scoring internal
mismatches, you are actually scoring a loop or bulge. To score terminal
mismatches, use `terminal_doublet`.

=item $o->terminal_doublet($t5, $t3, $b3, $b5, $left_side)

The arguments represent the same as above in addition to
C<left_side>, which is a boolean value that signifies if the
terminal end of the helix is on the left side.

This scores a doublet of RNA residues based on 
how well they would bind, assuming the doublet is found at a
terminal end of a helix structure. Then, it returns the score.
Terminal pairs must be at the end of both strands.

=item $o->bind($seq, $rRNA)

force_bind() forces an alignment between the two strands, 
assuming that they both pair together from their first bases on.  
It ignores both loops and bulges as possibilities.  Only bases 
in this forced alignment that can form helices are scored, and only the
score of the best helix is returned. That is, if there are two
or more helices formed, separated by gaps, only the score of the
lowest scoring helix is returned.

This method returns the total amount of free energy given off by
allowing two sequences of nucleotides to bind to each other.

For example, C<scan_brightly> prints the return values for all
subsequences found in a gene sequence and bound to the rRNA sequence
to produce the free energy calculations needed for the GSPtools model.

=back

=head1 COPYRIGHT AND LICENSE

This code is part of Kidnap. Kidnap originates from Joshua
Starmer's free2bind project as of July 2007. Original code
copyright Joshua Starmer, 2004. All modifications copyright
Hao Lian, 2007.

Kidnap is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Help me, I'm trapped Richard Stallman's hair.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 SEE ALSO

These papers contain the constants needed to implement the parameter
modules.

=over

=item Freier et al., PNAS (1986) Vol. 83, pp. 9373-9377

=item Jaeger et al., PNAS (1989) Vol. 86, pp. 7706-7710

=item SantaLucia, John, PNAS (1998) Vol. 95, pp. 1460-1465

=item Xia et al., Biochemistry (1998), Vol. 37, pp. 14719-14735

=item Mathews et al., J. Mol. Biol. (1999), Vol. 288, pp. 911-940

=back

=cut