use warnings; use strict;
use Smooth qw(prot2codon codon2prot);
package Dysentery;

# Takes an array (not a reference to an array) of
# codons--the array @_--and prints them neatly in
# rows of ten.
sub print_run {
    for (my $i = 0; $i < @_; $i += 10) {
        my $end = $i + 9;
        # In the common case the last codon doesn't
        # end on a multiple of 10, Perl will slice
        # some undefs for us.
        $end = @_ % 10 + $i - 1 if $end >= @_;
        
        my @codons = @_[$i .. $end];
        print join ' ', @codons, $/;
    }
}

# Returns a fabulous list of randomly selected codons
# that jive with the given protein sequence. Automates
# what the frameshift kids did on the board on Thursday
# and Friday.
sub run {
    if ($#ARGV == 1)
        { die "Perhaps you meant to use --check.$/" }
    my @codons = ();
    Smooth::webopen shift @ARGV, sub {
        # Sanitize and check for empty lines, upon which
        # I return.
        local $_ = $_;
        s/[^A-Z]//g;
        return unless $_;
        # Parse each and every character and store.
        map { push @codons, prot2codon($_) } split //;
    };
    # Arbitrary stop codon and 12-leader sequence
    # chosen by Vivek.
    print 'gcc aua ggc uau', $/;
    push @codons, 'uga'; @codons;
}

# Returns 1 if the codon sequence matches the
# protein sequence. Else returns 0.
#
# In addition, it returns both protein
# sequences.
sub check {
    ##### Get a list of codons.
    my @codons = (); my $proteins;
    Smooth::webopen shift @ARGV,
        sub { push @codons, split /\s+/ };
    
    # -2 because the last codon is the stop codon.
    @codons = @codons[4 .. @codons-2];
    map { $proteins .= codon2prot($_) } @codons;
    
    ##### Get actual list of proteins.
    my $actual = '';
    Smooth::webget shift @ARGV,
        sub { $actual = join '', @_ };
    $actual =~ s/[^A-Z]//g;
    
    return (1,$proteins,$actual) if $actual eq $proteins;
    (0,$proteins,$actual);
}

# Takes a truth value and prints the news accordingly.
sub print_check {
    my ($truth, $given, $actual) = @_;
    ##### Output results with strings for eye comparison
    ##### since wdiff is out of the question.
    if ($truth) { print "$/# They are EQUAL. Congratulations.$/" }
    else { print "$/# They are NOT EQUAL.$/" }    
    print " Given: $given$/$/";
    print "Actual: $actual$/";
}

if ($0 eq __FILE__) {
    my $help = <<END;
NAME
    dysentery.pl

USAGE
    1) Convert protein sequence (1-letter abbreviation standard) to codons.
        `perl dysentery.pl [url/path to protein sequence]`
    2) Check if a codon sequence produces a protein sequence.
        `perl dysentery.pl [url/path to codons] [url/path to proteins]`
END
    if (!@ARGV or $ARGV[0] eq '--help') { print $help; }
    elsif ($ARGV[0] eq '--check') { shift @ARGV; print_check check; }
    else { print_run run; }
}
1;