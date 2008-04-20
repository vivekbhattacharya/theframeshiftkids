use strict; use warnings;

open(my $h, shift);
my @lines = <$h>;
close($h);

sub infer {
    my ($i, $start, $stop) = @_;

    # Lucky for us, the next line is always the /gene. Then we'll
    # grope around in the dark and try to find /product.
    $i++ until $lines[$i] =~ /gene=/;
    my ($gene) = $lines[$i] =~ /gene="(.*?)"/;

    # Some products span more than one line. We'll need to remove
    # initial whitespace and newlines before extracting $desc.
    $i++ until $lines[$i] =~ /product=/;
    my $desc;
    for (join '::', @lines[$i .. $i+3]) {
        s/::\s+//g; s/(\r|\n)/ /g;
        ($desc) = /product="(.*?)"/;
    }
    return "$gene\t$start\t$stop\t$desc";
}

# Doing a foreach loop over <$h> (linear) and parsing (decidedly
# nonlinear) is hard, so we'll sacrifice elegance for sanity.
foreach my $i (0 .. $#lines) {
    local $_ = $lines[$i];
    next unless /^\s+CDS/;
    my ($start, $stop) = /(\d+)\.\.(\d+)/;

    # We're matching for complement(#..#). It's good to know, so
    # we'll mark it with a negative $start.
    $start *= -1 if /complement\(/;
    print infer($i, $start, $stop), "\n";
}
