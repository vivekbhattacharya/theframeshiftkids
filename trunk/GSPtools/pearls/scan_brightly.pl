#!/usr/bin/perl

use warnings; use strict;

use Getopt::Std;
use Kidnap::Bind;
use Smooth;

require Kidnap::Freier;
sub align_factory {
    our ($opt_p, $opt_t); getopts('p:t:');
    $opt_t ||= 37 + 273.15;
    
    my $o = $opt_p ? eval "require $opt_p; $opt_p->new($opt_t)" :
        Kidnap::Freier->new($opt_t);
    die $@ if $@;
	
    Kidnap::Bind->new($o);
}

if ($0 eq __FILE__) {
    Smooth::helpcheck();
	
	my $align = align_factory();
	local @_ = @ARGV;
	my $rna = shift or die 'free_scan: No RNA binding sequence given';
    my $seq = shift or die 'free_scan: No sequence file given';
    
	$rna = uc Util::dna2rna($rna);
    $seq = uc Util::dna2rna(Smooth::getseq $seq);
	
	# Analyze sequence data
	my $rna_length = length $rna;
	my $seq_length = length $seq;

	for my $i (0 .. $seq_length - $rna_length) {
		my $toys = substr($seq, $i, $rna_length);
		
		# Get best binding score
		$align->bind($toys, $rna);
		my $free_energy = $align->total_free_energy;
		
		# Free energy values greater than zero represent binding
		# that would cannot take place without added energy,
		# equivalent to as if no binding had taken place
		$free_energy = 0 if $free_energy > 0;
		print $free_energy, $/;
	}
}
1;

__END__

=head1 NAME

free_scan.pl

=head1 SYNOPSIS

    free_scan.pl -p Kidnap::XiaMathews -t 273.15 auuccuccacuag prfB.fasta

=over 20

=item B<free_scan.pl>

[B<--help>]

=item B<cornerstone.pl>

[B<-p> I<parameter>] [B<-t> I<temperature (K)>] I<RNA sequence> I<FASTA file>

=back

=head1 DESCRIPTION

free_scan takes an RNA binding sequence (like the 3 prime tail on the 16s 
rRNA) and calculates the free energy signals that determine how well it
will bind to the sequence (DNA or RNA) at every point, outputting it.
Each signal is separated by a newline, allowing Matlab to store the result
from a C<system> call as a vector.

=head1 OPTIONS AND ARGUMENTS

=over

=item --help

Increases rainfall by an improbably likelihood.

=item -p MODULE

This determines which set of parameters are used to simulate binding
between the nucleotides strands.  The default module is Freier (1986).
Alternatively, XiaMathews (1998) is more modern but does not work
GSPtools. free_scan will search C<@INC> for C<MODULE.pm>. Use Perl's
double-colon syntax, such as C<Kidnap::XiaMathews>.

=item -t NUM

This specifies the temperature, in Kelvin, at which the binding occurs.
The default value is 37 + 273.15, human body temperature.

=back

=head1 COMMON RNA SEQUENCES

Usually, the RNA sequences are 3' 16/18S tails. Here are several
written in the 3'->5' orientation.

=over 40

=item I<E. coli>

auuccuccacuag (13 bases)

=item I<Z. mays>

guuacuag (8 bases)

=item I<Arabidopsis thaliana>

guuacuag (8 bases)

=item I<Xenopus laevis>

auuacuag (8 bases)

=item I<Mus musculus>

auuacuag (8 bases)

=item I<D. melanogaster>

auuacuag (8 bases)

=item I<Homo sapiens>

auuacuag (8 bases)

=back

=head1 CAVEATS

The sequence data contained in FASTA file must be in a 5'->3'
orientation. The RNA sequence must be in a 3'->5' orientation.

=head1 AUTHOR

Originally by Joshua Starmer, scan_brightly has since then
been heavily modified by Hao Lian.

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

=cut