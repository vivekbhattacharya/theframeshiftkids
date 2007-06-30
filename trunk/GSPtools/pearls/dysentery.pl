# For testing:
# Codon sequence (rpoS): http://shadytrees.pastebin.ca/raw/594753
# Alternative: http://shadytrees.pastebin.ca/raw/594746
# Protein sequence: http://shadytrees.pastebin.ca/raw/595330

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
        s/[^A-Z]//g;
        return unless $_;
        map { push @codons, prot2codon($_) } split //;
    };
    # Arbitrary stop codon and 12-leader sequence
    use Data::Dumper;
    print $ARGV[1], $/;
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
    my $seq = getseq shift @_;
    # Split the sequence into threes using a nasty
    # trick involving exclamation points and intrigue.
    $seq =~ s/(...)/$1!/g;
    my @codons = map { codon2prot $_ }
        split('!', $seq);
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

our $help = <<END;
NAME
    dysentery.pl
    (web-enabled)
    
USAGE
    1) Convert protein sequence (1-letter abbreviation standard) to codons.
        `perl dysentery.pl [path to protein sequence] "[leader sequence]"`
    2) Check if a codon sequence produces a protein sequence.
        `perl dysentery.pl --check [path to codons] [path to proteins]`
    3) Check if two codon sequences produce the same protein sequence.
        `perl dysentery.pl --rcheck [path to codons] [path to codons]`
END
if ($0 eq __FILE__) {
    if (!@ARGV or $ARGV[0] eq '--help') { print $help; }
    # Shift out the --argument so it's business as usual.
    elsif ($ARGV[0] eq '--check') {
        print_check check($ARGV[1], $ARGV[2]);
    }
    elsif ($ARGV[0] eq '--rcheck') {
        print_check rcheck($ARGV[1], $ARGV[2]);
    }
    else { print_run run($ARGV[0]); }
}
1;