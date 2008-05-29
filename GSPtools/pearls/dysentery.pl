# For testing:
use warnings; use strict;
package Dysentery;
use Smooth qw(prot2codon codon2prot getseq);

# Takes an array (not a reference to an array) of
# codons--the array @_--and prints them neatly in
# rows of ten.
sub print_run {
    my $i = 0;
    foreach (@_) {
        print $_[$i], ' ';
        print $/ if ($i+1) % 10 == 0;
        $i++;
    }
}

# Takes a truth value and prints the news accordingly.
sub print_check {
    my ($truth, $given, $actual) = @_;
    # Output results with strings for eye comparison
    # since wdiff is out of the question.
    if ($truth) { print "$/# They are equal. Congratulations.$/" }
    else { print "$/# They are not equal.$/" }    
    print " Given: $given$/$/";
    print "Actual: $actual$/";
}

# Returns a fabulous list of randomly selected codons
# that jive with the given protein sequence. Automates
# what the frameshift kids did on the board on Thursday
# and Friday.
sub run {
    my @codons = ();
    Smooth::webopen shift, sub {
        # Sanitize and check for empty lines.
        # Then parse every character and store.
        return if /^>/;
        s/[^A-Z]//g;
        return unless $_;
        push @codons, prot2codon(split //);
    };
    # Arbitrary stop codon and 12-leader sequence
    use Data::Dumper;
    print $ARGV[2], $/;
    push @codons, 'uga'; @codons;
}

# Returns 1 if the codon sequence matches the
# protein sequence. Else returns 0.
#
# In addition, it returns both protein
# sequences.
sub check {
    # Remove 12-leader sequence and stop codons.
    my ($seq, $actual) = @_;
    my $proteins = substr(seq2proteins($seq), 4);
    $proteins =~ s/\.//g;

    $actual = Smooth::webslurp($actual);
    $actual =~ s/[^A-Z]//g;
    
    return (1,$proteins,$actual) if $actual eq $proteins;
    (0,$proteins,$actual);
}

# Returns a list of proteins, converted from
# a gene sequence.
sub seq2proteins {
    my @codons = codon2prot Smooth::seq2codons getseq shift;
    return join('', @codons);
};

# Converts two codon sequences into their respective
# protein sequences and returns a truth value.
sub rcheck {
    my ($left, $right) =
        (seq2proteins(shift), seq2proteins(shift));
    return (1,$left,$right) if $left eq $right;
    return (0,$left,$right);
}

if ($0 eq __FILE__) {
    my %table = (
                 '--help' => \&Smooth::helpcheck,
                 '--check' => sub { print_check check(@ARGV[1..2]) },
                 '--rcheck' => sub { print_check rcheck(@ARGV[1..2]) },
                 '--run' => sub { print_run run($ARGV[1]) }
                 );
    if (@ARGV > 0 and defined $table{$ARGV[0]}) {
        $table{$ARGV[0]}->();
    }
    else {
        Smooth::helpcheck();
    }
}
1;

__END__

=head1 NAME

dysentery.pl (web-enabled)

=head1 SYNOPSIS

dysentery.pl --run protein-sequence.txt "leader"

dysentery.pl --check codons.txt protein-sequence.txt

dysentery.pl --rcheck codons1.txt codons2.txt

=over 20

=item B<--run>

Converts a protein sequence using the standard 1-letter abbreviations
to codons

=item B<--check>

Checks if a codon sequence matches the protein sequence

=item B<--rcheck>

Checks if two codon sequences produces the same protein sequences

=back

=head1 URLs

These will help you test this program.

=over 20

=item rpoS codon sequence

http://shadytrees.pastebin.ca/raw/594753

=item Alternative rpoS codon sequence

http://shadytrees.pastebin.ca/raw/594746

=item Protein sequence

http://shadytrees.pastebin.ca/raw/595330

=back

=cut
