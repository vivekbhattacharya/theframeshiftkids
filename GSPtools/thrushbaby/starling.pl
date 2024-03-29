use strict; use warnings;
package Starling;
use File::Basename;
use lib dirname(__FILE__) . '/../pearls';
use Smooth qw(prot2codon codon2prot getseq);

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
    my ($loc, $file, $times) = @_;
    my $seq = getseq($file);
    
    # Fuzzy search (-3) to account for +1/-1 frameshifts.
    # Also, ignore the 12-leader (+12) sequence.
    # $i tells us how far we underestimated the location.
    #
    # Also, pull back (modulo) to the nearest *actual*
    # codon, ignoring all shifting.
    my $break = $loc*3 + 12;
    return strsmash($seq, $break, $times*3);
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
use List::Util qw(shuffle);
use File::Path qw(mkpath);
if ($0 eq __FILE__) {
    Smooth::helpcheck;
    
    my ($file, $folder, $loc) = @ARGV;
    # Increase 4 to a higher value if you're feeling
    # bold and you have a (very) fast computer.
    my ($before, $critical, $after) = beforehand $loc, $file, 4;
    
    my $i = 1;
    $folder .= "/$loc";
    mkpath($folder);
    map {
        open(my $handle, ">$folder/$i.txt");
        print $handle ($before . $_ . $after);
        $i++;
    } shuffle @{seq2permutations $critical};
}
1;

__END__

=head1 NAME

starling.pl (web-enabled)

=head1 SYNOPSIS

=over 20

=item B<starling.pl>

I<offset> I<sequence file> I<destination> I<codon> [I<codon> ...]

=back

=head2 EXAMPLE

    starling.pl "J:\chill\0\rpoS.txt" "J:\temp" 322
    
    starling.pl "J:\chill\179\98.txt" "J:\temp" 350

=head1 

=head1 OPTIONS AND ARGUMENTS

=over

=item I<offset>

starling.pl will search from this point onward in the sequence file.
If you're using the index starling.pl outputs, just feed that back
into this argument for the next one.

=item I<sequence file>

A list of codons that may or may not be in a FASTA format.

=item I<destination>

starling.pl, if the folder does not exist, creates it. In addition,
into it starling.pl outputs all the codon permutations under the format
of C<$folder/$location>. For example, "J:\temp\25" if the targeted
codon is the 25th codon.

=item I<codon> [I<codon> ...]

starling.pl parses this non-empty list of codons and locations.
These three metrics occur in the format C<$codon_$location>.
starling.pl will select the farthest neighbor of the earliest codon
in the list that occurs after the offset.

A codon Y is the neighbor of codon X if Y is a neighbor of a neighbor
of X or if Y is one codon away from X. For example, given the list
"uuu_25 ggg_26 aaa_28," starling.pl will--aside from creating the
permutations--choose "ggg_26" as the farthest neighbor and output
C<{26 'ggg,26'}> for Matlab to C<eval>uate.

There is one more caveat to this process, described below.

=back

=head1 CAVEATS

starling.pl chooses the earliest codon in the +0 reading frame
to correct any frameshifts that may have occurred during the
Matlab process, thus incorrectly reporting the list of codons. Thus,
dysentery.pl should correctly report any output file of starling.pl
equivalent to the original file inputted into thrushbaby.m.

=cut