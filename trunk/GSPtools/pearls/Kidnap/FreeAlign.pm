use warnings; use strict;
package Kidnap::FreeAlign;
## This module does free energy calculations.
##
## Authors: Joshua Starmer <jdstarme@unity.ncsu.edu>, 2004
##	The Frameshift Kids, 2007
## Copyright (C) 2004, The Frameshift Kids
##
## This product is licensed under the GNU General Public License
## <http://www.gnu.org/copyleft/gpl.html>.
##
## Values for the various parameters are from:
## 	Freier et al., PNAS (1986) Vol. 83, pp. 9373-9377
## 	Jaeger et al., PNAS (1989) Vol. 86, pp. 7706-7710
## 	SantaLucia, John, PNAS (1998) Vol. 95, pp. 1460-1465
## 	Xia et al., Biochemistry (1998), Vol. 37, pp. 14719-14735
## 	Mathews et al., J. Mol. Biol. (1999), Vol. 288, pp. 911-940

# Define global constants
use constant TRUE => 1;
use constant FALSE => 0;
use constant FAIL => -1;

use constant BIG_NUM => 30000;

# Binding and alignment
use constant FREIER => 1; # RNA parameters (1986)
use constant XIA_MATHEWS => 3; # RNA parameters (1998)

sub new {
	my $self = {};
    my ($class, $cat) = @_;
	# Are terminal/internal bulges allowed?
	$self->{NoBulge} = FALSE;
	$self->{InternalBulge} = FALSE;
	
	# Allow loops?
	$self->{Loop} = FALSE;
	
	# Binding temperature (K), related to dH and dS
	$self->{Temp} = 37 + 237.15;
	
	# Best score for binding two sequences together
	$self->{BestScore} = 0;
    $self->{Cat} = $cat;
	
	# Lengths for trace_route
	$self->{Seq} = {x => 0, y => 0};
	bless($self, shift);
}

## Returns the total amount of free energy given off by
## allowing two sequences of nucleotides to bind to each other.
sub total_free_energy {
	my $self = shift;
    $self->{BestScore} + $self->{Cat}->{InitPenalty};
}

## Returns free energy score of a "doublet" of RNA residue based on 
## how well they would bind, assuming that the doublet is found inside a
## helix structure. Thus, this only considers Watson/Crick paris plus G/U
## pairs.   A "doublet" consists of two pairs of RNA residues.
##
## NOTE: This only scores matches or internal (to a helix structure, i.e.
## not in a loop or bulge) G/U mismatches.  If you're scoring internal
## mismatches, you are actually scoring a loop or bulge. To score terminal
## mismatches, use `terminal_doublet`.
##
## t5: top strand, 5'
## t3: top strand, 3'
## b3: bottom strand, 3'
## b5: bottom strand, 5'
sub internal_doublet {
	my $self = shift;
    $self->{Cat}->internal_doublet(@_);
}


## terminal_doublet() scores a doublet[1] of RNA residues based on 
## how well they would bind, assuming the doublet is found at a
## terminal end of a helix structure. Then, it returns the score.
## Terminal pairs must be at the end of both strands.
## 	[1]: A doublet is two pairs of RNA residues.
## 
## ARGUMENTS: (the usual and)
## left_side: Is the terminal end of the helix on the left side?
sub terminal_doublet {
	my $self = shift;
    return $self->{Cat}->terminal_doublet(@_);
}

## Defines some substr calls used in force_bind. See
## force_bind() for more details.
##
## Note: fivethree is threefive for a y doublet and
## similarly ->{context} is switched.
package Strand;
use Data::Dumper;
# Asterisk needed to workaround substr's handling
# of negative numbers. That is, we want $i-1 to
# still work even when it doesn't produce s55.
sub new {
	my ($class, $seq, $length) = @_;
	my $self = {seq => '*'.$seq, length => $length};
	
    #helper($self, $i);
	bless $self, $class;
}

sub update {
    my ($self, $i) = @_;
    my ($s55, $s5, $s3, $s33) = split '', substr($self->{seq}, $i, 3);
	$self->{fivethree} = [$s5, $s3];
	$self->{context} = ($i > 0 && $i < $self->{length}-2) ? 
		[$s55, $s33] : ['', ''];
}

sub all { @{shift->{fivethree}} }
sub context { @{shift->{context}} }

## force_bind() forces an alignment between the two strands, 
## assuming that they both pair together from their first bases on.  
## It ignores both loops and bulges as possibilities.  Only bases 
## in this forced alignment that can form helices are scored, and only the
## score of the best helix is returned. That is, if there are two
## or more helices formed, separated by gaps, only the score of the
## lowest scoring helix is returned. The score does not include
## the one-time penalty for initiating a helix. (Use InitPenalty.)
##
## We also return the the position in the sequences where the best sub-helix
## begins (best_start) and its length (helix_length).
package Kidnap::FreeAlign;
sub bind {
    my ($self, $seq_x, $seq_y) = @_;    
    $self->{BestScore} = 0;

    my $length = length $seq_x;
	# I need at least two bases in each sequence in order to form
	# a structure between them.
    if ($length != length $seq_y) {
		die 'force_bind(): Unequal sequence lengths';
    } if ($length < 2) {
		die 'force_bind(): Sequences are too short';
    }

    my $best_score = my $score = 0;
    my $helix_start = my $helix_end = 0;
    my $best_start = -1;
    
    my $x = Strand->new($seq_x, $length);
    my $y = Strand->new($seq_y, $length);
    for my $i (0 .. $length-2) {
        $x->update($i); $y->update($i);
		
		# Start each helix structure with a "terminal doublet."
		if ($score == 0) {
			$score = $self->terminal_doublet($x->all, $y->all, 1);
		} else {
			# Some internal "doublets" need more context to
			# be scored correctly.  See the scoring for GU
			# in doublet_xia_mathews() for more details.	
			$score += $self->internal_doublet($x->all, $y->all, $x->context, $y->context);
		}
		
		if ($score > 0) {
			$helix_start = $i + 1;
			$score = 0;
		} elsif ($score < $best_score) {
			$best_start = $helix_start;
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
		
		my $internal_score = $self->internal_doublet($x->all, $y->all);
		my $terminal_score = $self->terminal_doublet($x->all, $y->all, !1);

		$best_score -= $internal_score - $terminal_score;
	
		# Here is where we should check for "self symmetry." With
		# the force_bind, we are assuming that the strands are
		# stretched out, and thus, not stuck to themselves forming
		# their own hairpin loop.  My understanding is that the "self symmetry"
		# penalty results from the individual strands being stuck together 
		# in their own hairpins.
    }

    $self->{BestScore} = $best_score;
    my $helix_length = ($best_start == -1) ?
        ($helix_end - $best_start + 2) : 0;
    ($best_score, $best_start, $helix_length);
}
1;