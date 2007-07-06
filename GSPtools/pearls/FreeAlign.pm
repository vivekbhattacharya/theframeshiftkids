use warnings; use strict;
package FreeAlign;
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
	# Are terminal/internal bulges allowed?
	$self->{NoBulge} = FALSE;
	$self->{InternalBulge} = TRUE;
	
	# Allow loops?
	$self->{Loop} = TRUE;
	
	# Binding temperature (K), related to dH and dS
	$self->{Temp} = 37 + 237.15;
	
	# Best score for binding two sequences together
	$self->{BestScore} = 0;
	
	# Lengths for trace_route
	$self->{Seq} = {x => 0, y => 0};
	bless($self, shift);
}

## Specify which set of binding parameters that you wish to use.
## Pass one of the constants above as the identifier.
sub select_parameters {
	my $self = shift;
    $self->{Parameter} = shift;
	
	# This penalty is one-time for forming any bonds between
	# two different RNA strands (Lewin, Genes VI, page 102).
	# See Freier (1986).
    if ($self->{Parameter} == FREIER) { $self->{InitPenalty} = 3.4; }
	
	# Helix initiation penalty, see Xia (1998).
	elsif ($self->{Parameter} == XIA_MATHEWS) {
		my ($dH, $dS) = (3.61, -1.5);
		$self->{InitPenalty} = $dH -
			$self->{Temp}*($dS/1000);
    }
}

## Returns the total amount of free energy given off by
## allowing two sequences of nucleotides to bind to each other.
sub total_free_energy {
	my $self = shift;
    if (($self->{Parameter} == FREIER) ||
		($self->{Parameter} == XIA_MATHEWS)) {
		return ($self->{BestScore} + $self->{InitPenalty});
    }
}

## Returns DNA or RNA complement.
sub complement {
    my ($self, $seq, $is_rna) = @_;
    
    if ($is_rna) { $seq =~ tr/aucgAUCG/uagcUAGC/; }
	else { $seq =~ tr/atcgATCG/tagcTAGC/; }
	return $seq;
}


sub dna2rna {
    my ($self, $seq) = @_;
    $seq =~ tr/tT/uU/;
	return $seq;
}

## Calculates and returns the free energy (score) for an internal loop
## given the number of nucleotides in the loop. This is similar
## to a mismatch penalty in a normal alignment.
sub internal_loop {
    my ($self, $length) = @_;

    # are we allowed to consider loops in an optimal binding?    
    if (!$self->{Loop}) { return BIG_NUM; }
    if ($length == 0) { return 0; }
    if ($length == 1) {
		# Single base bulge
		print STDERR "internal_loop: Loop length of one is not a real loop";
		return BIG_NUM;
    }

	for ($self->{Parameter}) {
		if ($_ == FREIER || $_ == XIA_MATHEWS) {
			return $self->loop_jaeger($length);
		}
	}
}

sub loop_jaeger {
    my ($self, $l) = @_;

	# Short loops (<10 bases) are scored independent of temperature.
	# See Jaeger (1989), where these are experimentally derived.
	my %scores = (2 => 4.1, 3 => 4.5, 4 => 4.9, 5 => 5.1,
		6 => 5.7, 7 => 5.9, 8 => 6.0, 9 => 6.1, 10 => 6.3);
	if ($l >= 2 && $l <= 10) { return $scores{$l} }
	
	return 5.7 + $self->jaeger_formula($l/6);
}

# Jaegar 1989
sub jaeger_formula {
	my ($self, $length) = @_;
	return (1.75 * 1.987/1000 * $self->{Temp} * log($length));
}

## Calculates and returns the free energy for a bulge given
## the length of nucleotides in the loop.
## This is similar to a gap penalty in a normal alignment.
##
## When scoring by hand, remember that the bases in the bulge
## are ignored and you calculate an "doublet" for the surrounding
## bases.
sub bulge_energy {
    my ($self, $l) = @_;
    # Are we allowed to consider bulges in an optimal binding?
    if (!$self->{InternalBulge}) { return BIG_NUM; }
	
	my %scores = (0 => 0, 1 => 1, 2 => 3.1, 3 => 3.5,
		4 => 4.2, 5 => 4.8, 6 => 5.0, 7 => 5.2, 8 => 5.3,
		9 => 5.4, 10 => 5.5);
    if ($l >= 0 && $l <= 10) { return $scores{$l}; }
	return 4.8 + $self->jaeger_formula($l/5 * 310.15); 
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
    for ($self->{Parameter}) {
		return $self->doublet_freier(@_) if $_ == FREIER;
		return $self->doublet_xia_mathews(@_) if $_ == XIA_MATHEWS;
    }
}

# Freier 1986
sub doublet_freier {
    my ($self, $t5, $t3, $b3, $b5) = @_;
    if (!$self->valid_pair($t5, $b3) ||
		!$self->valid_pair($t3, $b5)) { return BIG_NUM; }

	my %scores = (
		'AU' => {
			# Watson/Crick matches
			'AU' => -0.9, 'UA' => -0.9, 'GC' => -1.7, 'CG' => -2.1,
			# G/U mismatches
			'GU' => -0.5, 'UG' => -0.7
		}, 'UA' => {
			'AU' => -1.1, 'UA' => -0.9, 'GC' => -1.8, 'CG' => -2.3,
			'GU' => -0.7, 'UG' => -0.5
		}, 'CG' => {
			'AU' => -1.8, 'UA' => -1.7, 'GC' => -2.9, 'CG' => -2.0,
			'GU' => -1.5, 'UG' => -1.5
		}, 'GC' => {
			'AU' => -2.3, 'UA' => -2.1, 'GC' => -2.9, 'CG' => -3.4,
			'GU' => -1.3, 'UG' => -1.9
		}, 'GU' => {
			'AU' => -0.5, 'UA' => -0.7, 'GC' => -1.5, 'CG' => -1.9,
			'GU' => -0.5, 'UG' => -0.5
		}, 'UG' => {
			'AU' => -0.7, 'UA' => -0.5, 'GC' => -1.5, 'CG' => -1.3,
			'GU' => -0.6, 'UG' => -0.5
		},
	);
    my $score = $scores{$t5 . $b3}->{$t3 . $b5};
    return $score if $score;
    return BIG_NUM;
}


## SantaLucia (1998)
sub doublet_santalucia {
    my ($self, $t5, $t3, $b3, $b5) = @_;
    # We are scoring DNA, therefore only score WC pairs.
    if (!$self->wc_pair($t5, $b3) ||
		!$self->wc_pair($t3, $b5)) {
		return BIG_NUM;
    }
	
	my %scores = (
		'AT' => {
			'AT' => -1.00, 'TA' => -0.88, 'GC' => -1.28, 'CG' => -1.44
		}, 'TA' => {
			'AT' => -0.58, 'TA' => -1.00, 'GC' => -1.45, 'CG' => -1.30
		}, 'CG' => {
			'AT' => -1.45, 'TA' => -1.28, 'GC' => -2.17, 'CG' => -1.42
		}, 'GC' => {
			'AT' => -1.30, 'TA' => -1.44, 'GC' => -1.42, 'CG' => -1.84
		},
	);
    my $score = $scores{$t5 . $b3}->{$t3 . $b5};
    return $score if $score;
    return BIG_NUM;
}


## Xia (1998) (Watson-Crick pairs) and Mathews (1999) (G/U)
sub doublet_xia_mathews {
	my $self = shift;
    my ($t5, $t3, $b3, $b5, $t55, $t33, $b33, $b55) = @_;
    if (!$self->valid_pair($t5, $b3) ||
		!$self->valid_pair($t3, $b5)) {
		return BIG_NUM;
    }
    
	my %scores = (
		'AU' => {
			'AU' => [-6.82, -19.0], 'UA' => [-9.38, -26.7],
			'GC' => [-10.48, -27.1], 'CG' => [-11.40, -29.5],
			'GU' => [-3.21, -8.6], 'UG' => [-8.81, -24.0],
		}, 'UA' => {
			'AU' => [-7.69, -20.5], 'UA' => [-6.82, -19.0],
			'GC' => [-10.44, -26.9], 'CG' => [-12.44, -32.5],
			'GU' => [-6.99, -19.3], 'UG' => [-12.83, -37.3],
		}, 'CG' => {
			'AU' => [-10.44, -26.9], 'UA' => [-10.48, -27.1],
			'GC' => [-10.64, -26.7], 'CG' => [-13.39, -32.7],
			'GU' => [-5.61, -13.5], 'UG' => [-12.11, -32.2],
		}, 'GC' => {
			'AU' => [-12.44, -32.5], 'UA' => [-11.40, -29.5],
			'GC' => [-13.39, -32.7], 'CG' => [-14.88, -36.9],
			'GU' => [-8.33, -21.9], 'UG' => [-12.59, -32.5],
		}, 'GU' => {
			'GU' => [-13.47, -41.82], 'UG' => [-14.59, -51.2],
			'UA' => [-8.81, -24.0], 'AU' => [-12.83, -37.3],
			'GC' => [-12.11, -32.3], 'CG' => [-12.59, -32.5],
		}, 'UG' => {
			'GU' => [-9.26, -30.8], 'UG' => [-13.47, -41.82],
			'UA' => [-3.21, -8.6], 'AU' => [-6.99, -19.3],
			'GC' => [-5.61, -13.5], 'CG' => [-8.33, -21.9],
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



## terminal_doublet() scores a doublet[1] of RNA residues based on 
## how well they would bind, assuming the doublet is found at a
## terminal end of a helix structure. Then, it returns the score.
## Terminal pairs must be at the end of both strands.
## 	[1]: A doublet is two pairs of RNA residues.
## 
##	ARGUMENTS:
## 		t5: top strand, 5'
## 		t3: top strand, 3'
## 		b3: bottom strand, 3'
## 		b5: bottom strand, 5'
##		left_side: Is the terminal end of the helix on the left side?
sub terminal_doublet {
	my $self = shift;
    for ($self->{Parameter}) {
		return $self->terminal_doublet_freier(@_) if $_ == FREIER;
		return $self->terminal_doublet_xia_mathews(@_) if $_ == XIA_MATHEWS;
    }
}

sub terminal_doublet_xia_mathews {
	my $self = shift;
    my $doublet_score = $self->doublet_xia_mathews(@_);

	my ($t5, $t3, $b3, $b5, $left_side) = @_;
	my @chesthair = $left_side ? ($t5, $b3) : ($t3, $b5);
	my $terminal_penalty = $self->terminal_pair_xia_mathews(@chesthair);

    return $doublet_score + $terminal_penalty;
}


## This scores a pair of RNA residues based on theoretical binding
## ability, assuming the pair exists at the terminal end of a helix
## structure. Terminal pairs must be at the ends of both strands.
sub terminal_pair_xia_mathews {
    my ($self, $t, $b) = @_;

    my ($dH, $dS) = (0, 0);
	if (grep $t.$b, ('AU', 'UA', 'GU', 'UG')) {
		($dH, $dS) = (3.72, 10.5);
    }
    return $dH - $self->{Temp}*($dS/1000);
}

sub terminal_doublet_freier {
	# Terminal pairs must be at the very ends of both strands.
    my ($self, $t5, $t3, $b3, $b5, $left_side) = @_;

	# Flip bases if we are looking at the left side of the helix.
    if ($left_side) {
		($t5, $b5) = ($b5, $t5);
		($t3, $b3) = ($b3, $t3);
    }

    # We only need the first pair to match, provided it
	# is not a G/U mismatch.
    if (!$self->wc_pair($t5, $b3)) { return BIG_NUM; }
	my %scores = (
		'AU' => {
			# Watson/Crick matches
			'AU' => -0.9, 'UA' => -0.9, 'GC' => -1.7, 'CG' => -2.1,
			# G/U mismatches
			'GU' => -0.9, 'UG' => -0.9,
			# Mismatches
			'AA' => -0.8, 'CC' => -0.7, 'GG' => -1.0, 'UU' => -0.8,
			'AC' => -1.0, 'AG' => -1.0, 'CA' => -0.7, 'GA' => -0.8,
			'UC' => -0.8, 'CU' => -0.7,
		}, 'UA' => {
			'AU' => -1.1, 'UA' => -0.9, 'GC' => -1.8, 'CG' => -2.3,
			'GU' => -0.9, 'UG' => -1.0,
			'AA' => -1.0, 'CC' => -0.6, 'GG' => -1.2, 'UU' => -0.5,
			'AC' => -0.8, 'AG' => -1.1, 'CA' => -0.7, 'GA' => -1.1,
			'UC' => -0.6, 'CU' => -0.5,
		}, 'CG' => {
			'AU' => -1.8, 'UA' => -1.7, 'GC' => -2.0, 'CG' => -2.9,
			'GU' => -1.6, 'UG' => -1.9,
			'AA' => -1.9, 'CC' => -1.1, 'GG' => -1.9, 'UU' => -1.2,
			'AC' => -2.0, 'AG' => -1.9, 'CA' => -1.0, 'GA' => -1.0,
			'UC' => -1.5, 'CU' => -0.8,
		}, 'GC' => {
			'AU' => -2.3, 'UA' => -2.1, 'GC' => -2.9, 'CG' => -3.4,
			'GU' => -1.4, 'UG' => -2.3,
			'AA' => -1.1, 'CC' => -0.6, 'GG' => -1.4, 'UU' => -0.7,
			'AC' => -1.7, 'AG' => -1.3, 'CA' => -1.1, 'GA' => -1.6,
			'UC' => -0.8, 'CU' => -0.5,
		},
	);
	my $score = $scores{$t5 . $b3}->{$t3 . $b5};
	return $score if $score;
	return BIG_NUM;
}


## This scores 3' dangling RNA residues at the ends of helices
## based on theoretical binding ability. The score is for the
## unpaired terminal nucleotide.
##
## ARGUMENTS:
##	t5: top strand, 5'
##	t3: top strand, 3'
##	b3: bottom strand, paired with t5.
sub dangling_3pdrime {
    my ($self, $t5, $t3, $b3) = @_;

    if ($self->{Parameter} == FREIER) {
		return $self->dangling_3prime_freier($t5, $t3, $b3);
    } elsif ($self->{Parameter} == XIA_MATHEWS) {
		# CURRENTLY USING THE FREIER PARAMETERS... SHOULD USE PARAMETERS
		# FOUND IN...
		# Serra and Turner, 1995 - Predicting thermodyamic properties of RNA.
		#  Methods in Enzymology 259, 242-261
		return $self->dangling_3prime_freier($t5, $t3, $b3);
    }
}

sub dangling_3prime_freier {
    my ($self, $t5, $t3, $b3) = @_;

    # We need the Watson-Crick pair to match.
    return BIG_NUM if !$self->wc_pair($t5, $b3);
    # Do we allow terminal bulges?
    return BIG_NUM if $self->{NoBulge};
	
	my %scores = (
		'AU' => {
			'A' => -0.8, 'C' => -0.5, 'G' => -0.8, 'U' => -0.6,
		}, 'CG' => {
			'A' => -1.7, 'C' => -0.8, 'G' => -1.7, 'U' => -1.2,
		}, 'GC' => {
			'A' => -1.1, 'C' => -0.4, 'G' => -1.3, 'U' => -0.6,
		}, 'UA' => {
			'A' => -0.7, 'C' => -0.1, 'G' => -0.7, 'U' => -0.1,
		},
	);
	my $score = $scores{$t5 . $b3}->{$t3};
	return $score if $score;
    return BIG_NUM;
}


sub dangling_5prime {
    my ($self, $t5, $t3, $b5) = @_;

    # we need the pair to match (and no G/U mismatches)...
    if (!($self->wc_pair($t3, $b5))) {
	return BIG_NUM;
    }    
    
    # are we allowing terminal bulges?
    if ($self->{NoBulge}) {
	return BIG_NUM;
    }
    
    if (($t3 eq 'A') && ($b5 eq 'U')) {
	if ($t5 eq 'A') {                # AA
	    return -0.3; # Freier (1986)    U
	}
	if ($t5 eq 'C') {                # CA
	    return -0.3; # Freier (1986)    U
	}
	if ($t5 eq 'G') {                # GA
	    return -0.4; # Freier (1986)    U
	}
	if ($t5 eq 'U') {                # UA
	    return -0.2; # Freier (1986)    U
	}
    }

    if (($t3 eq 'C') && ($b5 eq 'G')) {
	if ($t5 eq 'A') {                # AC
	    return -0.5; # Freier (1986)    G
	}
	if ($t5 eq 'C') {                # CC
	    return -0.2; # Freier (1986)    G
	}
	if ($t5 eq 'G') {                # GC
	    return -0.2; # Freier (1986)    G
	}
	if ($t5 eq 'U') {                # UC
	    return -0.1; # Freier (1986)    G
	}
    }

    if (($t3 eq 'G') && ($b5 eq 'C')) {
	if ($t5 eq 'A') {                # AG
	    return -0.2; # Freier (1986)    C
	}
	if ($t5 eq 'C') {                # CG
	    return -0.3; # Freier (1986)    C
	}
	if ($t5 eq 'G') {                # GG
	    return -0.0; # Freier (1986)    C
	}
	if ($t5 eq 'U') {                # UG
	    return -0.0; # Freier (1986)    C
	}
    }

    if (($t3 eq 'U') && ($b5 eq 'A')) {
	if ($t5 eq 'A') {                # AU
	    return -0.3; # Freier (1986)    A
	}
	if ($t5 eq 'C') {                # CU
	    return -0.2; # Freier (1986)    A
	}
	if ($t5 eq 'G') {                # GU
	    return -0.2; # Freier (1986)    A
	}
	if ($t5 eq 'U') {                # UU
	    return -0.2; # Freier (1986)    A
	}
    }

    print STDERR "ERROR: FreeAlign::dangling_5prime()\n";
    print STDERR "ERROR: 5' dangle $t5$t3\n";
    print STDERR "                    $b5 not scored yet...\n";	    
    return BIG_NUM;
}


## Returns true if two characters form a valid helix pair,
## i.e. a standard Watson-Crick pair or a G-U mismatch.
sub valid_pair {
    my ($self, $a, $b) = @_;
	my $ab = $a . $b;
	return scalar grep /$ab/, qw/GC CG AU AT UA TA GU UG/;
}


## Returns true if two characters form a valid W-C pair
sub wc_pair {
    my ($self, $a, $b) = @_;
	my $ab = $a . $b;
	return scalar grep /$ab/, qw/GC CG AU AT UA TA/;
}

## Defines some substr calls used in force_bind. See
## force_bind() for more details.
package Strand;
sub new {
	my ($class, $seq, $i, $length) = @_;
	my $self = {};
	
	$self->{fivethree} = [substr($seq, $i, 1), substr($seq, $i+1, 1)];
	$self->{context} =
		($i < $length-2) ? 
		[substr($seq, $i-1, 1), substr($seq, $i+2, 1)] :
		['', ''];
	
	bless $self, $class;
}

sub all {
	my $self = shift;
	return @{$self->{fivethree}};
}

sub context {
	my $self = shift;
	return @{$self->{context}};
}

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
package FreeAlign;
sub force_bind {
    my ($self, $seq_x, $seq_y) = @_;    
    $self->{BestScore} = 0;

    my $length = length($seq_x);
	# I need at least two bases in each sequence in order to form
	# a structure between them.
    if ($length != length($seq_y)) {
		die 'force_bind(): Unequal sequence lengths', $/;
    } if ($length < 2) {
		die 'force_bind(): Sequences are too short', $/;
    }

    my $score = 0;
    my $best_score = 0;

    my $helix_start = 0;
    my $best_start;
    my $helix_end = 0;
    for my $i (0 .. $length-2) {
		my $x = Strand->new($seq_x, $i, $length);
		my $y = Strand->new($seq_y, $i, $length);
		
		# Start each helix structure with a "terminal doublet."
		if ($score == 0) {
			$score = $self->terminal_doublet($x->all, $y->all, TRUE);
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
			$helix_end = $i;
			$best_score = $score;
			$best_start = $helix_start;
		}
    }

	# Swap out the last internal for a terminal doublet score
	# so that the helix starts and ends with a terminal doublet.
	# We have to do this now because we cannot know in advance
	# where this terminal doublet is going to be.
    if ($best_score < 0) {
		my $x = Strand->new($seq_x, $helix_end, 0, 0);
		my $y = Strand->new($seq_y, $helix_end, 0, 0);
		
		my $internal_score = $self->internal_doublet($x->all, $y->all);
		my $terminal_score = $self->terminal_doublet($x->all, $y->all, FALSE);

		$best_score -= $internal_score - $terminal_score;
	
		# Here is where we should check for "self symmetry." With
		# the force_bind, we are assuming that the strands are
		# stretched out, and thus, not stuck to themselves forming
		# their own hairpin loop.  My understanding is that the "self symmetry"
		# penalty results from the individual strands being stuck together 
		# in their own hairpins.
    }

    $self->{BestScore} = $best_score;
    my $helix_length =
		(defined $best_start) ? ($helix_end - $best_start + 2) : 0;
    return ($best_score, $best_start, $helix_length);
}
1;