use strict; use warnings;
package Starling;
use Smooth qw(prot2codon codon2prot getseq);

sub parse {
    return sort { $a->[2] <=> $b->[2] }
        map {
            my ($codon, $loc, $times) = split '_';
            [$times, $codon, $loc];
        } @_;
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

sub pick {
    my $start = shift;
    printf q/{%s '%s,%s' %s}/, $_[0]->[2], $_[0]->[1];
    $_[0];
}

sub seq2permutations {
    my $seq = shift;
    my @proteins = codon2prot Smooth::seq2codons $seq;
    Smooth::permute map {
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
    of the four[2] codons behind the most "common"
    frameshifter after the start, including of
    course all the codons surrounding that area.
    
    [1]: In the form of codon_location_occurrence
    [2]: Anything more than four is expensive but
        possible Contact the authors if you have
        an extremely fast machine :).
END
if ($0 eq __FILE__) {
    if (!@ARGV or $ARGV[0] eq '--help') { print $help; exit }
    
    my ($start, $file, $folder, @rest) = @ARGV;
    my @data = parse @rest;
    my ($critical, $before, $after) = beforehand pick($start, @data), $file, 2;
    
    my $i = 1;
    map {
        open(my $handle, ">$folder/$i.txt");
        print $handle ($before . $_ . $after);
        $i++;
    } @{seq2permutations $critical};
}
1;