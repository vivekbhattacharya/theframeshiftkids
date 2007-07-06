use strict; use warnings;
package Starling;
use Smooth qw(prot2codon codon2prot getseq);

sub parse {
    my $start = shift;
    return
        grep { $_->[2] > $start + 9 }
        sort { $a->[2] <=> $b->[2] }
        map {
            my ($codon, $loc, $times) = split '_';
            [$times, $codon, $loc];
        } @_;
}

# Parses and returns the string segments
# before and after the needle in the
# following diagram:
#       xxxxxxx-->xxxxxxxx<--xxxxx
#
# $str := Haystack
#   $a := Start of the needle
#   $b := Length of the needle
# Returns (before, needle, after)
sub strsmash {
    my ($str, $a, $b) = @_;
    return (
        substr($str, 0, $a-$b),
        substr($str, $a-$b, $b),
        substr($str, $a),
    );
}

# Return the codon sequence the comes before
# the given codon. For example, uga's
# sequence would be ...aggggguaucuu for prfB.
#
# $data: pointer to triplet
# $file: path to file
# $times: number of codons
sub beforehand {
    my ($data, $file, $times) = @_;
    my $seq = getseq($file);
    
    # Fuzzy search (-3) to account for +1/-1 frameshifts.
    # Also, ignore the 12-leader (+12) sequence.
    # $i tells us how far we underestimated the location.
    #
    # Also, pull back (modulo) to the nearest *actual*
    # codon, ignoring all shifting.
    my $break = $data->[2]*3 + 12 - 3;
    for (substr $seq, $break) {
        $break += index $_, $data->[1];
        $break -= $break % 3;
        return strsmash($seq, $break, $times*3);
    }
}

sub pick {
    # Time to end ThrushBaby.
    if ($#_ < 0) { print q/{-1 -1}/; exit }
    
    # Gobble up neighboring frameshifts
    # if they're only one away.
    my $data = shift;
    foreach my $wombat (@_) {
        if ($wombat->[2] - $data->[2] == 1) {
            $data = $wombat;
        } else { +last }
    }
    
    # Communicate with Matlab ... now!
    my (undef, $codon, $loc) = @$data;
    printf q/{%s '%s,%s'}/, $loc, $codon, $loc;
    return $data;
}

sub permute {
    use List::Util qw(reduce);
    no warnings 'once';
    # Note: we have to store a temporary
    # $x to access $_ after it's been shadowed.
    reduce {[
        map {
            my $x = $_;
            map "$_$x", @$a;
        } @$b
    ]} @_;
}

sub seq2permutations {
    my $seq = shift;
    my @proteins = codon2prot Smooth::seq2codons $seq;
    permute map {
        [split /,/, $Smooth::expression{$_}]
    } @proteins;
}
        
use Data::Dumper;
our $help = <<END;
NAME
    jovial.pl
    (part of ThrushBaby)

USAGE
    jovial.pl (1) (2) (3) (4)
        1) Start position for `pick`ing
        2) Gene sequence file
        3) Work folder (already existing)
        4) and a list of triplets[1]
    
    I will fill the work folder with permutations
    of the four[2] codons behind[3] the earliest[4]
    frameshifter after the start, including of
    course all the codons surrounding that area.
    
    [1]: In the form of codon_location_occurrence
    [2]: Anything more than four is expensive but
        possible contact the authors if you have
        an extremely fast machine :).
    [3]: Where behind is behind after correcting
        any possible shifts that may have occurred
        so that the replacement codons do not actually
        change yielded protein.
    [4]: Technically, jovial considers neighboring
        codons as one single group and chooses the
        farthest right neighbor.[5]
    [5]: A neighbor to a codon is either one away
        from the codon or one away from another
        neighbor.
END
if ($0 eq __FILE__) {
    if (!@ARGV or $ARGV[0] eq '--help') { print $help; exit }
    
    my ($start, $file, $folder, @rest) = @ARGV;
    my @data = parse $start, @rest;
    # Increase 4 to a higher value if you're feeling
    # bold and you have a (very) fast computer.
    my ($before, $critical, $after) = beforehand pick(@data), $file, 4;
    
    my $i = 1;
    map {
        open(my $handle, ">$folder/$i.txt");
        print $handle ($before . $_ . $after);
        $i++;
    } @{seq2permutations $critical};
}
1;