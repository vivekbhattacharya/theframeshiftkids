#!/usr/bin/perl

##========================================================================
##
## Author: Joshua Starmer <jdstarme@unity.ncsu.edu>, 2004
## 
## Copyright (C) 2004, Joshua Starmer
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
## 
##========================================================================

use warnings; use strict;

use Getopt::Std; 
use FreeAlign;

my $sHelp = <<END_OF_HELP;

Usage: free_scan.pl [options] rna_seq sequence.fasta
    
    free_scan takes an RNA binding sequence (like the 3 prime tail on the 16s 
    rRNA) and sees how well it will bind to the sequence (DNA or RNA, if DNA
    it is first converted to RNA) at each point in a range specified in 
    index_data.
    
    All of the sequence contained in sequence.fasta will be examined.

    NOTES:
    The sequence data contained in "sequence.fasta" is assumed to be
    in 5'->3' orientation.

    "rna_seq" is assumed to be written in 3'->5' orientation.

    example:
    shell> ./free_scan.pl -p XIA_MATHEWS auuacuagg human_dna.fasta human_gene_indices > delta_G_values.txt

    Commonly used strings for "rna_seq" are 3 prime 16/18 S tails.
    Here are several listed as reference, written 3prime->5prime:

    E. coli: auuccuccacuag (13 bases)
    Z. mays: guuacuag (8 bases)
    Arabidopsis thaliana: guuacuag (8 bases)
    Xenopus laevis: auuacuag (8 bases)
    Mus musculus: auuacuag (8 bases)
    D. melanogaster: auuacuag (8 bases)
    H. sapiens: auuacuag (8 bases)

Options:

-h
   print this help message.

-p FREIER|XIA_MATHEWS
   Allows you to determine which set of parameters are used to simulate binding
   between the to strands of nucleotides.  The default value is XIA_MATEHWS.
   FREIER: RNA binding parameters from 1986
   XIA_MATHEWS: RNA binding parameters from 1998

-t temperature
   The temperature, in Kelvin, at which the binding is supposed to take place.
   The default value is 37 + 273.15.
END_OF_HELP


use Smooth;
if ($0 eq __FILE__) {
    Smooth::helpcheck();

    my %opts;
    getopts('ehdqcL:T:n:fp:t:Ors:g', \%opts);
     
    my $align = FreeAlign->new();
	$opts{'p'} ||= 'FREIER';

	for ($opts{'p'}) {
		$align->select_parameters($align->FREIER), last if $_ eq 'FREIER';
		$align->select_parameters($align->XIA_MATHEWS), last if $_ eq 'XIA_MATHEWS';
		die 'free_scan: Invalid doublet ID parameter (-p)';
	}

    $align->{Temp} = $opts{'t'} if defined $opts{'t'};

	local @_ = @ARGV;
	my $rna = shift or die 'free_scan: No RNA binding sequence given';
    my $seq = shift or die 'free_scan: No sequence file given';
    
	$rna = uc $align->dna2rna($rna);
    $seq = uc $align->dna2rna(Smooth::getseq $seq);
	
	# Analyze sequence data
	my $rna_length = length $rna;
	my $seq_length = length $seq;

	for my $i (0 .. $seq_length - $rna_length) {
		my $seq = substr($seq, $i, $rna_length);
		
		# Get best binding score
		$align->force_bind($seq, $rna);
		my $free_energy = $align->total_free_energy;
		
		# Free energy values greater than zero represent binding
		# that would cannot take place without added energy,
		# equivalent to as if no binding had taken place
		$free_energy = 0 if $free_energy > 0;
		print $free_energy, $/;
	}
}
1;