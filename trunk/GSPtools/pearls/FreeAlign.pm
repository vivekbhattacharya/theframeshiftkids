use warnings; use strict;
package FreeAlign;

##========================================================================
## This is a perl module or class that implements some common 
## subroutines that free energy calculations by dynamic programming may want
## to use
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
## See later comments
##========================================================================

# DEFINE GLOBAL CONSTANTS
use constant TRUE => 1;
use constant FALSE => 0;
use constant FAIL => -1;

use constant DEBUG => FALSE;
use constant BIG_NUM => 30000;

# CONSTANTS FOR PRETTY PRINTING VARIOUS THINGS
use constant MATRIX_CELL_WIDTH => 8; # good for floating point numbers...
use constant LINE_LENGTH => 77; # good for xterms

# CONSTANTS FOR SCORING THE BINDING/ALIGNMENT:
use constant FREIER => 1; # RNA parameters (1986)
use constant SANTALUCIA => 2; # DNA (NOT RNA) parameters (1998)
use constant XIA_MATHEWS => 3; # RNA parameters (1998)

# HELIX_START_GC and HELIX_START_AT are both for DNA hybridization
# and should be added on for both ends.  If the helix starts with 
# and ends with a G/C pair, then I should use HELIX_START_GC twice.
# See Santa Lucia (1998).
use constant HELIX_START_GC => 0.98; 
use constant HELIX_START_AT => 1.03;

# SINGLE_STRANDED_END corresponds to a terminal gap penalty in alignment
# terminology.
use constant SINGLE_STRANDED_END => 0;
use constant MIN_LOOP_STEP => 2;

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
		} elsif ($_ == SANTALUCIA) {
			return $self->loop_santalucia($length);
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


## Returns free enegrey score of a "doublet" of RNA residue based on 
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
		return $self->doublet_santalucia(@_) if $_ == SANTALUCIA;
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
    ## We are scoring DNA, therefore only score WC pairs.
    if (!($self->wc_pair($t5, $b3) && $self->wc_pair($t3, $b5))) {
	return BIG_NUM;
    }
    
    if (($t5 eq "A") && ($b3 eq "T")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "T")) { # AA
	    return -1;    # SantaLucia (1998) TT
	}
	if (($t3 eq "T") && ($b5 eq "A")) { # AT
	    return -0.88; # SantaLucia (1998) TA
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # AG
	    return -1.28; # SantaLucia (1998) TC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # AC
	    return -1.44; # SantaLucia (1998) TG
	}
    }
    
    if (($t5 eq "T") && ($b3 eq "A")) {
	# watson/crick matches...
	if (($t3 eq "T") && ($b5 eq "A")) { # TT
	    return -1;    # SantaLucia (1998) AA
	}
	if (($t3 eq "A") && ($b5 eq "T")) { # TA
	    return -0.58; # SantaLucia (1998) AT
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # TG
	    return -1.45; # SantaLucia (1998) AC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # TC
	    return -1.30; # SantaLucia (1998) AG
	}
    }
    
    if (($t5 eq "C") && ($b3 eq "G")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "T")) { # CA
	    return -1.45; # SantaLucia (1998) GT
	}
	if (($t3 eq "T") && ($b5 eq "A")) { # CT
	    return -1.28; # SantaLucia (1998) GA
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # CC
	    return -1.42; # SantaLucia (1998) GG
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # CG
	    return -2.17; # SantaLucia (1998) GC
	}
    }
    
    if (($t5 eq "G") && ($b3 eq "C")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "T")) { # GA
	    return -1.30; # SantaLucia (1998) CT
	}
	if (($t3 eq "T") && ($b5 eq "A")) { # GT
	    return -1.44; # SantaLucia (1998) CA
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # GG
	    return -1.42; # SantaLucia (1998) CC
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # GC
	    return -1.84; # SantaLucia (1998) CG
	}       
    }
    
    print STDERR "ERROR: FreeAlign::doublet_santalucia()\n";
    print STDERR "ERROR: doublet $t5$t3\n";
    print STDERR "               $b3$b5 not scored yet...\n";
    return BIG_NUM;
    
}


##========================================================================
##
## SUBROUTINE doublet_xia_matthews(), parameter values are from
##   Xia (1998) (for Watson-Crick pairs) and Mathews (1999) (for G/U 
##   mismatches).
##
##   NOTE: These parameter values are generated on the fly...
##
##   See "doublet" for more details.
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub doublet_xia_mathews {
    my ($self, $t5, $t3, $b3, $b5, $t55, $t33, $b33, $b55) = @_;
    
    if (!($self->valid_pair($t5, $b3) && $self->valid_pair($t3, $b5))) {
	return BIG_NUM;
    }
    
    my $deltaH;
    my $deltaS;
    
    if (($t5 eq "A") && ($b3 eq "U")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # AA
	    $deltaH = -6.82;  # Xia (1998)    UU
	    $deltaS = -19.0;  # Xia (1998)
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # AU
	    $deltaH = -9.38;  # Xia (1998)    UA
	    $deltaS = -26.7;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # AG
	    $deltaH = -10.48; # Xia (1998)    UC
	    $deltaS = -27.1;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # AC
	    $deltaH = -11.40; # Xia (1998)    UG
	    $deltaS = -29.5;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  AG
	    $deltaH = -3.21; # Mathews (1999)  UU
	    $deltaS = -8.6;  # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  AU
	    $deltaH = -8.81; # Mathews (1999)  UG
	    $deltaS = -24.0; # Mathews (1999)
	}
    }
    
    elsif (($t5 eq "U") && ($b3 eq "A")) {
	# watson/crick matches...
	if (($t3 eq "U") && ($b5 eq "A")) { # UU
	    $deltaH = -6.82;  # Xia (1998)    AA
	    $deltaS = -19.0;  # Xia (1998)
	}
	if (($t3 eq "A") && ($b5 eq "U")) { # UA
	    $deltaH = -7.69;  # Xia (1998)    AU
	    $deltaS = -20.5;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # UG
	    $deltaH = -10.44; # Xia (1998)    AC
	    $deltaS = -26.9;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # UC
	    $deltaH = -12.44; # Xia (1998)    AG
	    $deltaS = -32.5;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  UG
	    $deltaH = -6.99; # Mathews (1999)  AU
	    $deltaS = -19.3; # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  UU
	    $deltaH = -12.83; # Mathews (1999) AG 
	    $deltaS = -37.3;  # Mathews (1999)
	}
    }
    
    elsif (($t5 eq "C") && ($b3 eq "G")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # CA
	    $deltaH = -10.44; # Xia (1998)    GU
	    $deltaS = -26.9;  # Xia (1998)
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # CU
	    $deltaH = -10.48; # Xia (1998)    GA
	    $deltaS = -27.1;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # CC
	    $deltaH = -13.39; # Xia (1998)    GG
	    $deltaS = -32.7;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # CG
	    $deltaH = -10.64; # Xia (1998)    GC
	    $deltaS = -26.7;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  CG
	    $deltaH = -5.61; # Mathews (1999)  GU
	    $deltaS = -13.5; # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  CU
	    $deltaH = -12.11; # Mathews (1999) GG
	    $deltaS = -32.2;  # Mathews (1999)
	}
    }
    
    elsif (($t5 eq "G") && ($b3 eq "C")) {
	# watson/crick matches...
	if (($t3 eq "A") && ($b5 eq "U")) { # GA
	    $deltaH = -12.44; # Xia (1998)    CU
	    $deltaS = -32.5;  # Xia (1998)
	}
	if (($t3 eq "U") && ($b5 eq "A")) { # GU
	    $deltaH = -11.40; # Xia (1998)    CA
	    $deltaS = -29.5;  # Xia (1998)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { # GG
	    $deltaH = -13.39; # Xia (1998)    CC
	    $deltaS = -32.7;  # Xia (1998)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { # GC
	    $deltaH = -14.88; # Xia (1998)    CG
	    $deltaS = -36.9;  # Xia (1998)
	}
	
	# G/U mismatches...
	if (($t3 eq "G") && ($b5 eq "U")) { #  GG
	    $deltaH = -8.33; # Mathews (1999)  CU
	    $deltaS = -21.9; # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  GU
	    $deltaH = -12.59; # Mathews (1999) CG
	    $deltaS = -32.5;  # Mathews (1999)
	}
    } 

    elsif (($t5 eq "G") && ($b3 eq "U")) {	
	if (($t3 eq "G") && ($b5 eq "U")) { #  GG
	    $deltaH = -13.47; # Mathews (1999) UU
	    $deltaS = -41.82; # Based on Mathews (1999) - see note in Table 4
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  GU
	    $deltaH = -14.59; # Mathews (1999) UG
	    $deltaS = -51.2;  # Mathews (1999)
	    
	    if (defined($t55) && defined($t33) &&     # only scoring middle two
		defined($b33) && defined($b55)) {     #  /\
		if (($t55 eq "G") && ($b33 eq "C") && # GGUC
		    ($t33 eq "C") && ($b55 eq "G")) { # CUGG
		    
		    # Mathews (1999) scores the entire GGUC/CUGG unit with
		    # deltaH = -30.80 and deltaS = -86.0.  The values below
		    # are calculated by... 
		    # deltaH = -30.80-(2*-8.33), -8.33 is the deltaH for GG/CU
		    # deltaS = -86-(2*-21.9), -21.9 is the deltaS for GG/CU
		    
		    $deltaH = -14.14; # Based on Mathews (1999)
		    $deltaS = -42.2;  # Based on Mathews (1999)
		}
	    }	
	}
	if (($t3 eq "U") && ($b5 eq "A")) { #  GU
	    $deltaH = -8.81;  # Mathews (1999) UA
	    $deltaS = -24.0;  # Mathews (1999)	
	}
	if (($t3 eq "A") && ($b5 eq "U")) { #  GA
	    $deltaH = -12.83; # Mathews (1999) UU
	    $deltaS = -37.3;  # Mathews (1999)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { #  GG
	    $deltaH = -12.11; # Mathews (1999) UC
	    $deltaS = -32.2;  # Mathews (1999)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { #  GC
	    $deltaH = -12.59; # Mathews (1999) UG
	    $deltaS = -32.5;  # Mathews (1999)
	}
    }    

    elsif (($t5 eq "U") && ($b3 eq "G")) {	
	if (($t3 eq "G") && ($b5 eq "U")) { #  UG
	    $deltaH = -9.26;  # Mathews (1999) GU
	    $deltaS = -30.8;  # Mathews (1999)
	}
	if (($t3 eq "U") && ($b5 eq "G")) { #  UU
	    $deltaH = -13.47; # Mathews (1999) UU
	    $deltaS = -41.82; # Based on Mathews (1999) - see note in Table 4
	}
	if (($t3 eq "U") && ($b5 eq "A")) { #  UU
	    $deltaH = -3.21;  # Mathews (1999) GA
	    $deltaS = -8.6;   # Mathews (1999)
	}
	if (($t3 eq "A") && ($b5 eq "U")) { #  UA
	    $deltaH = -6.99;  # Mathews (1999) GU
	    $deltaS = -19.3;  # Mathews (1999)
	}
	if (($t3 eq "G") && ($b5 eq "C")) { #  UG
	    $deltaH = -5.61;  # Mathews (1999) GC
	    $deltaS = -13.5;  # Mathews (1999)
	}
	if (($t3 eq "C") && ($b5 eq "G")) { #  UC
	    $deltaH = -8.33;  # Mathews (1999) GG
	    $deltaS = -21.9;  # Mathews (1999)
	}
    }

    if (!defined($deltaH) || !defined($deltaS)) {
	print STDERR "ERROR: FreeAlign::doublet_xia_mathews()\n";
	print STDERR "ERROR: doublet $t5$t3\n";
	print STDERR "               $b3$b5 not scored yet...\n";
	return BIG_NUM;
    }

#    print STDERR "deltaH: $deltaH\n";
#    print STDERR "deltaS: $deltaS\n";

    return $deltaH - $self->{Temp}*($deltaS/1000);
}



##========================================================================
##
## SUBROUTINE terminal_doublet() scores a "doublet" of RNA residues based on 
##   how well they would bind, assuming that the doublet is found at a
##   terminal end of a helix structure and returns the score.
##   A "doublet" consists of two pairs of RNA residues.
##
##   NOTE: TERMINAL PAIRS MUST BE AT THE VERY ENDS OF BOTH STRANDS
##
##   EXAMPLES: GC would get one score and AU would get another...
##             CG                         UG
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   left_side - is the terminal end of the helix on the left side,
##     or on the right side?
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub terminal_doublet {
    my ($self, $t5, $t3, $b3, $b5, $left_side) = @_;

    if ($self->{Parameter} == FREIER) {
	return $self->terminal_doublet_freier($t5, $t3, $b3, $b5, $left_side);
    } elsif ($self->{Parameter} == SANTALUCIA) {
	return $self->terminal_doublet_santalucia($t5, $t3, 
						  $b3, $b5, $left_side);
    } elsif ($self->{Parameter} == XIA_MATHEWS) {
	return $self->terminal_doublet_xia_mathews($t5, $t3, 
						   $b3, $b5, $left_side);
    } else {
	print STDERR "ERROR: FreeAlign::terminal_doublet()\n";
	print STDERR "       gParameter: $self->{Parameter}, undefined\n";
	return BIG_NUM;
    }
}

##========================================================================
##
## SUBROUTINE terminal_doublet_xia_mathews(), parameter values are from
##   Xia (1998) and Mathews (1999)
##
##   See "terminal_dobulet" for more detials.
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, 3'
##
##   b5 - bottom strand, 5'
##
##   left_side - is the terminal end of the helix on the left side,
##     or on the right side?
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub terminal_doublet_xia_mathews {
    my ($self, $t5, $t3, $b3, $b5, $left_side) = @_;
    
    my $doublet_score = $self->doublet_xia_mathews($t5, $t3, $b3, $b5);

    my $terminal_penalty;
    if ($left_side) {
	$terminal_penalty = $self->terminal_pair_xia_mathews($t5, $b3);
    } else {
	$terminal_penalty = $self->terminal_pair_xia_mathews($t3, $b5);
    }

    return $doublet_score + $terminal_penalty;
}


##========================================================================
##
## SUBROUTINE terminal_pair_xia_mathews() scores a pair of RNA residues based 
##   on how well they would bind, assuming that the pair is found at a
##   terminal end of a helix structure and returns the score.
##
##   NOTE: TERMINAL PAIRS MUST BE AT THE VERY ENDS OF BOTH STRANDS
##
##   EXAMPLES: G would get one score and G would get another...
##             C                         U
##
##   ARGUMENTS:
##   t - top strand
##     
##   b - bottom strand
##
##   RETURN VALUES:
##   - a free energy value (a score) for the doublet.
##
##========================================================================
sub terminal_pair_xia_mathews {
    my ($self, $t, $b) = @_;

    my $deltaH = 0;
    my $deltaS = 0;

    if ((($t eq "A") && ($b eq "U")) ||
	(($t eq "U") && ($b eq "A")) ||
	(($t eq "G") && ($b eq "U")) ||
	(($t eq "U") && ($b eq "G"))) {
	$deltaH = 3.72; # Xia (1998) and Mathews (1999)
	$deltaS = 10.5; # Xia (1998) and Mathews (1999) 
    }
    
    return $deltaH - $self->{Temp}*($deltaS/1000);

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


##========================================================================
##
## SUBROUTINE dangling_3prime() scores 3' dangling RNA residues at the ends of
##   helices based on how well they would bind and returns the score.
##
##   EXAMPLES: GC would get one score and AU would get another...
##             C                          U
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b3 - bottom strand, paired with t5.
##
##   RETURN VALUES:
##   - a free energy value (a score) for the unpaired terminal nucleotide
##
##========================================================================
sub dangling_3prime {
    my ($self, $t5, $t3, $b3) = @_;

    if ($self->{Parameter} == FREIER) {
	return $self->dangling_3prime_freier($t5, $t3, $b3);
    } elsif ($self->{Parameter} == SANTALUCIA) {
	return $self->dangling_3prime_santalucia($t5, $t3, $b3);
    } elsif ($self->{Parameter} == XIA_MATHEWS) {
	# CURRENTLY USING THE FREIER PARAMETERS... SHOULD USE PARAMETERS
	# FOUND IN...
	# Serra and Turner, 1995 - Predicting thermodyamic properties of RNA.
	#  Methods in Enzymology 259, 242-261
	return $self->dangling_3prime_freier($t5, $t3, $b3);
    } else {
	print STDERR "ERROR: FreeAlign::dangling_3prime()\n";
	print STDERR "       gParameter: $self->{Parameter}, undefined\n";
	return BIG_NUM;
    }
}

## Returns score for a 3' dangling RNA residues at the 
## ends of helices based on how well they would bind.
## The score is for the unpaired terminal nucleotide.
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


##========================================================================
##
## SUBROUTINE dangling_5prime() scores 5' dangling RNA residues at the ends of
##   helices based on how well they would bind and returns the score.
##
##   EXAMPLES: GC would get one score and AU would get another...
##              C                          U
##
##   ARGUMENTS:
##   t5 - top strand, 5'
##     
##   t3 - top strand, 3'
##
##   b5 - bottom strand, paired with t3.
##
##   RETURN VALUES:
##   - a free energy value (a score) for the unpaired terminal nucleotide.
##
##========================================================================
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


##========================================================================
##
## SUBROUTINE force_bind() forces an alignment between the two strands, 
##   assuming that they both pair together from their first bases on.  
##   This method ignores both loops and bulges as possibilities.  Only bases 
##   in this forced alignment that can form helices are scored, and only the
##   score of the best helix is returned (that is to say, if there are two
##   or more helices formed, separated by gaps, then only the score of the
##   lowest scoring helix is returned.)
##
##   This is roughly
##   equivelent to assuming that both strands are pulled tight (imagine
##   holding two pieces of rope of equal length together by their ends and
##   pulling the ends apart from each other as far as you can).
##
##   NOTE: The sequences that are paired together are assumed to have the
##   same length.
##
##   ARGUMENTS:
##   seq_x - one of two sequences to be paired together.
##
##   seq_y - one of two sequences to be paired together.
##
##   RETURN VALUES:
##   best_score - The score for the lowest scoring helix in the binding.
##     For example, if the binding is...
##       GGGGCCCCCggAAUU
##       CCCCGGGGGggUUAA
##     Then only the score for the first helix (the one containing Gs and Cs)
##     is returned.
##     NOTE: This score DOES NOT INCLUDE the one time penalty for initiating
##     a helix.  This penalty can be subtracted later using the 
##     gInitPenalty value.
##
##   best_start - the position in the sequences where the best sub-helix
##     begins.
##
##   helix_length - the length of the best sub-helix.
##
##========================================================================
sub force_bind {
    my ($self, $seq_x, $seq_y) = @_;    

    $self->{BestScore} = 0;

    my $x_length = length($seq_x);
    my $y_length = length($seq_y);

    if ($x_length != $y_length) {
	print STDERR "FreeAlign::force_bind() - unequal sequence lengths.\n";
	return;
    }

    if (($x_length < 2) || ($y_length < 2)) {
	# you need at least two bases in each sequence in order to form
	# a structure between them.
	print STDERR "FreeAlign::force_bind() - sequences are too short.\n";
	return;	
    }


    my $current_score = 0;
    my $best_score = 0;

    my $helix_start_pos = 0;
    my $best_start;

    my $end_of_helix = 0;   

    for (my $i=0; $i<($x_length-1); $i++) {
	my $x_5 = substr($seq_x, $i, 1);
	my $x_3 = substr($seq_x, ($i+1), 1);
	
	my $y_3 = substr($seq_y, $i, 1);
	my $y_5 = substr($seq_y, ($i+1), 1);
	
	if ($current_score == 0) {
	    # start each helix structure with a "terminal doublet"
	    $current_score = $self->terminal_doublet($x_5, $x_3, 
						     $y_3, $y_5, TRUE);
	} else {
	    my $x_55 = "";
	    my $x_33 = "";
	    my $y_33 = "";
	    my $y_55 = "";

	    if ($i < $x_length-2) {
		# What is going on here is that some internal "doublets" need 
		# more context in order to be scored correctly.  See the 
		# scoring for GU in doublet_xia_mathews() for more details.
		#             UG
	        $x_55 = substr($seq_x, ($i-1), 1);
		$x_33 = substr($seq_x, ($i+2), 1);

		$y_33 = substr($seq_y, ($i-1), 1);
		$y_55 = substr($seq_y, ($i+2), 1);
	    }

	    $current_score = $current_score + 
		$self->internal_doublet($x_5, $x_3, 
					$y_3, $y_5, 
					$x_55, $x_33,
					$y_33, $y_55);
	}
	
	if ($current_score > 0) {
	    $helix_start_pos = $i + 1;
	    $current_score = 0;
	} elsif ($current_score < $best_score) {
	    $end_of_helix = $i;
	    $best_score = $current_score;
	    $best_start = $helix_start_pos;
	}
    }

    if ($best_score < 0) {
	# now swap out the last internal doublet score for a terminal doublet
	# score so that the helix starts and ends with a terminal doublet.
	# We have to do this now, since we can not know in advance
	# where this terminal doublet is going to be...

	my $x_5 = substr($seq_x, $end_of_helix, 1);
	my $x_3 = substr($seq_x, ($end_of_helix+1), 1);
	
	my $y_3 = substr($seq_y, $end_of_helix, 1);
	my $y_5 = substr($seq_y, ($end_of_helix+1), 1);
	
	my $internal_score = $self->internal_doublet($x_5, $x_3, $y_3, $y_5);
	
	my $terminal_score = $self->terminal_doublet($x_5, $x_3, 
						     $y_3, $y_5, FALSE);
    
	$best_score = $best_score - $internal_score + $terminal_score;

	# Here is where we should check for "self symmetry", however, given
	# the assumptions made with "force_bind", I'm not entirely sure it
	# is appropriate...  With the force_bind, we are assuming that the
	# strands are stretched out, and thus, not stuck to themselves forming
	# their own hairpin loop.  My understanding is that the "self symmetry"
	# penalty results from the individual strands being stuck together 
	# in their own hairpins.
    }

    $self->{BestScore} = $best_score;

    my $helix_length = 0;
    if (defined($best_start)) {
	$helix_length = $end_of_helix - $best_start + 2;
    }

    return ($best_score, $best_start, $helix_length);
}

##========================================================================
##
## SUBROUTINE print_matrix() prints out a matrix to STDERR
##
##   ARGUMENTS:
##   matrix - a reference to a matrix to print
##
##   width - the number of columns in the matrix
##
##   height - the number of rows in the matrix
##
##   row - a reference to an array of strings to label the rows
##
##   column - a reference to an array of strings to label the columns
##
##   floats - a boolean that indicates whether or not the matrix contains
##     floating point numbers or just integers.
##
##========================================================================
sub print_matrix {
    my ($self, $matrix, $width, $height, $row, $column, $floats) = @_;

    no warnings;

    my $print_string = "\%".MATRIX_CELL_WIDTH."d";

    if (defined($row)) {
	print STDERR " ";
	for (my $i=0; $i<$width; $i++) {
	    print STDERR " "x(MATRIX_CELL_WIDTH-1).$$row[$i];
	}
	print STDERR "\n";
    }

    for (my $y=0; $y<$height; $y++) {
	if (defined($column)) {
	    print STDERR $$column[$y];
	}

	for (my $x=0; $x<$width; $x++) {
	    if (defined($floats) && $floats) {
		my $temp = sprintf("\%.1f", $$matrix[$x][$y]);
		$print_string = " "x(MATRIX_CELL_WIDTH - length($temp)).$temp;
	    }
	    printf(STDERR $print_string, $$matrix[$x][$y]);
	}
	print STDERR "\n";
    }
}



##========================================================================
##
## SUBROUTINE print_binding():  pretty prints out the binding
##   between the two RNA strands.
##
##   ARGUMENTS:
##   seq_x - a sequence that has been bound to another sequence
##
##   x_label - the name of seq_x
##
##   seq_y - a sequence that has been bound to another sequence
##
##   y_label - the name of seq_y
##
##   match_string - an optional string that is intended to go between the
##     two sequences.  It is intended to show relationships between the two 
##     sequences (exact matches, close matches...)
##
##========================================================================
sub print_binding {
    my ($self, $seq_x, $x_label, $seq_y, $y_label, $match_string) = @_;

    my $label_length;
    my $x_label_length = length($x_label);
    my $y_label_length = length($y_label);
    if ($x_label_length > $y_label_length) {
	$label_length = $x_label_length;
	# pad $y_label with white space so that it is the same width as
	# $x_label...
	my $pad = $x_label_length - $y_label_length;
	$y_label = $y_label." "x$pad;
    } else {
	$label_length = $y_label_length;
	# pad $x_label with white space to make it the same width as 
	# $y_label...
	my $pad = $y_label_length - $x_label_length;
	$x_label = $x_label." "x$pad;	
    }

    my $alignment_length = length($seq_x);
    my $index_length = length($alignment_length);
    my $index_str = sprintf("%%%dd", $index_length);
    my $index_pad = " "x$index_length;

    my $line_length = LINE_LENGTH - $label_length - $index_length - 2;
    my $match_title = " "x$label_length;
    my $i=0;
    my $x_index = 0;
    my $y_index = 0;
    while ($i<$alignment_length) {
	my $x_sub = substr($seq_x, $i, $line_length);
	my $match_sub;
	if ($match_string) {
	    $match_sub = substr($match_string, $i, $line_length);
	}
	my $y_sub = substr($seq_y, $i, $line_length);	
	printf("%s $index_str: %s\n", $x_label, $x_index, $x_sub);
	if ($match_string) {	    
	    printf("%s %s: %s\n", $match_title, $index_pad, $match_sub);
	}
	printf("%s $index_str: %s\n\n", $y_label, $y_index, $y_sub);
	
	$x_index += ($x_sub =~ tr/a-zA-Z0-9//);
	$y_index += ($y_sub =~ tr/a-zA-Z0-9//);
	$i=$i+$line_length;
    }
}



##========================================================================
##
## SUBROUTINE read_fasta(): read a file in FASTA format.
##
##   NOTE: this function reads from wherever <$filehandle> last left off and
##     returns the next FASTA sequence.  To read all of the sequences in
##     a FASTA file, you may have to call this function multiple times.
##
##   ARGUMENTS:
##   filehandle - a reference to a filehandle to read the sequence from.  
##
##   RETURN VALUES:
##   seq - a sequence from the file.
##   
##   seq_accession - the accession number in the file.  NOTE: depending on
##     how the line of comments is formatted, this may or may not contain
##     the actual accession number.
##
##   seq_comments - the comments that come before the sequence.
##
##========================================================================
sub read_fasta {
    my ($self, $filehandle) = @_;

    # eliminate all of the filler until the first sequence
    my $input = "";
    my $done = FALSE;
    while ((!$done) && ($input = <$filehandle>)) {
	if ($input =~ /\A\x3e/) {  # this line starts with '>'
	    $done = TRUE;
	}
    }

    if (!$input) {
	# if there is nothing left in the file, we can just return 
	# everything undefined...
	return;
    }

    chomp($input);

    # anything on the same line as the first '>' character is a comment
    my $seq_comments = $input;
    # extract the accession # for the sequence from the first line of comments
    $input =~ /.*?\|.*?\|.*?\|(.*?)\|/;
    my $seq_accession = $1;

    # get the sequence
    my $seq;
    $done = FALSE;
    while((!$done) && ($input = <$filehandle>)) {
	if ($input =~ /\A\W/) { # this line starts with white space
	    $done = TRUE;
	} else {
	    chomp($input);
	    # eliminate any other whitespaces
	    $input =~ s/\s//;
	    $seq .= $input;
	}
    }
    return ($seq, $seq_accession, $seq_comments);
}


##========================================================================
##
## SUBROUTINE load_fasta(): opens and reads the first sequence in a FASTA file
##
##   NOTE: this function only returns the very first sequence in the FASTA
##     file.  If there are more than one sequence in the file and you want
##     to access them, you should pass the filehandle to "read_fasta()".
##
##   ARGUMENTS:
##   input_file - the name of the file to read the sequence from or it is a 
##     filehandle (like STDIN).  
##
##   RETURN VALUES:
##   seq - a sequence from the file.
##   
##   seq_accession - the accession number in the file.  NOTE: depending on
##     how the line of comments is formatted, this may or may not contain
##     the actual accession number.
##
##   seq_comments - the comments that come before the sequence.
##
##   filehandle - the filehandle for the newly opened file, or, if the 
##     input_file was a filehandle to begin with, input_file.
##
##========================================================================
sub load_fasta {
    my ($self, $input_file) = @_;
    
    my $input = "";
    my $done = FALSE;
    my $filehandle;
    if (-f $input_file) {
	open($filehandle, $input_file) || die("can't open $input_file: $!");
    } else {
	$filehandle = $input_file;
    }
    
    my ($seq, $seq_accession, $seq_comments) = $self->read_fasta($filehandle);
    
    return ($seq, $seq_accession, $seq_comments, $filehandle);
}


1;