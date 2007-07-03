use strict; use warnings;
package Jovial;
use Smooth qw(prot2codon codon2prot getseq);

sub parse {
    my $garbage = shift;
    return sort { $a->[2] <=> $b->[2] }
        map {
            my ($codon, $loc, $times) = split '_';
            [$times, $codon, $loc];
        } split /\s+/, $garbage;
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
    for (substr $seq, $breakpoint) {
        my $i = index $_, $data->[1];
        $breakpoint += $i;
        
        my $prefix = substr($seq, $breakpoint - $times*3, $times*3);
        my $before = substr($seq, 0, $breakpoint - $times*3);
        my $after = substr($seq, $breakpoint);
        return ($prefix, $before, $after);
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
    my ($file, $folder, $garbage) = @ARGV;
    my @data = parse $garbage;
    my ($critical, $before, $after) = beforehand $data[0], $file, 4;
    
    my $i = 1;
    map {
        open(my $handle, ">$folder/$i.txt");
        print $handle ($before . $_ . $after);
        $i++;
    } @{seq2permutations $critical};
}
1;