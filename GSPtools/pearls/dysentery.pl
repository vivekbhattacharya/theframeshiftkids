use Smooth;
use warnings;
use strict;

my @codons = ();
Smooth::webopen shift @ARGV, sub {
    my $line = $_;
    
    # http://en.wikipedia.org/wiki/List_of_standard_amino_acids#Gene_expression_and_biochemistry
    $line =~ s/[^A-Z]//g;
    my @chars = split //, $line;
    map { push @codons, Smooth::cupid($_); } @chars;
};
push @codons, 'uga';
print 'gcc aua ggc uau', $/;

my $i = 0;
my $total = scalar @codons;
while(1) {
    my $end = $i + 9;
    # In the common case the last codon doesn't
    # end on a multiple of 10, Perl will slice
    # some undefs for us.
    $end = $total % 10 + $i - 1 if $end > $total;
    
    my @subcodons = @codons[$i .. $end];
    print join ' ', @subcodons, $/;
    
    $i += 10;
    last if $i > @codons;
}