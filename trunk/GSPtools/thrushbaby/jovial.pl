use strict; use warnings;
package Jovial;
use Smooth qw(prot2codon codon2prot getseq);

sub parse {
    return sort { $a->[2] <=> $b->[2] }
        map {
            my ($codon, $loc, $times) = split '_';
            [$times, $codon, $loc];
        } @_;
}

# The Kissing Algorithm
# Imagine you have a bunch of codons with scores
# sitting single-file. Each codon receives the
# half the score of the following codon if the
# two codons are sitting next to each other in
# a gene sequence. This way, we account for
# very close frameshifts as the same problem area.
#
# Arguments: ([7, 'uga', 25], [3, 'ccc', 49], ...)
sub distribute_kisses {
    my $last = $_[0];
    # ->[0] is the count
    # ->[2] is the position
    for (my $i = 1; $i < @_; $i++) {
        my $now = $_[$i];
        if ($now->[2] - $last->[2] == 1) {
            my ($x, $y) = ($now->[0]/2, $last->[0]/2);
            $last->[0] += $x;
            $now->[0] += $y;
        }
        $last = $now;
    }
    @_;
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
    my $pos = $data->[2];
    my $seq = getseq($file);
    
    # Fuzzy search (-3) to account for +1/-1 frameshifts.
    # Also, ignore the 12-leader (+12) sequence.
    # $i tells us how far we underestimated the location.
    my $breakpoint = $pos*3 + 12 - 3;
    for (substr $seq, $breakpoint, -1) {
        my $i = index $_, $data->[1];
        my $prefix = substr($seq,
            $breakpoint + $i - $times*3, $times*3);
        return $prefix;
    }
}

sub seq2permutations {
    my $seq = shift;
    my @proteins = codon2prot Smooth::seq2codons $seq;
    Smooth::permute map {
        [split /,/, $Smooth::expression{$_}]
    } @proteins;
}
        
use Data::Dumper;
if ($0 eq __FILE__) {
    my @data = (['7','uga',25],['2','aac',26],['3','ccc',49],['5','uaa',57],['2','aug',146],['2','uga',171],['1','ucc',314],['1','cga',315]);
    my $preseq = beforehand $data[1], 'J:\frameshift\prfb.txt', 4;
    print Dumper(seq2permutations $preseq);
}
1;