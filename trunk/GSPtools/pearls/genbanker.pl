use strict; use warnings;
use ParseGenbank;
use Smooth;

use File::Path qw(mkpath);
use File::Spec;
use List::Util qw(max min);

my $genbank_file = shift;
my $genome_file = shift;
my $dir = shift or die 'No output directory given.';

my @lines = do {open(my $h, $genbank_file); <$h>};
sub infer {
    my ($i) = @_;

    $i++ until $lines[$i] =~ /gene=/;
    my ($gene) = $lines[$i] =~ /gene="(.*?)"/;
    my $gene_line = $i;

    # Genes can have a /product= xor a /pseudo line.
    my $pseudo = 0;
    for (; ; $i++) {
        last if $lines[$i] =~ /product=/;
        if ($lines[$i] =~ /\/pseudo/) {
            $pseudo = 1; last;
        }
    }

    my $desc;
    if ($pseudo) {$desc = '/pseudo';}
    # Some products span more than one line. We'll need to remove
    # initial whitespace and newlines before extracting $desc.
    else {
        for (join '::', @lines[$i .. $i+3]) {
            s/::\s+//g; s/(\r|\n)/ /g;
            ($desc) = /product="(.*?)"/;
        }
    }
    return ($gene_line, $gene, $desc);
}

my $genome = lc do {
    local $/;
    open(my $genome_h, $genome_file)
      or die 'Unable to open genome.';
    <$genome_h>;
};

$genome =~ s/^>.*$//m; # Get rid of FASTA header.
$genome =~ s/\s//g; # Get rid of newlines.
$genome =~ tr/t/u/; # Convert to mRNA.

# Write back a sanitized version for faster future runs.
open(my $genome_h2, ">$genome_file")
  or die 'Unable to open genome (for writing).';
print $genome_h2 $genome;
close($genome_h2);

unless (-d $dir) {
    mkpath $dir or die 'Unable to write to output directory.';
}

my $parser = new ParseGenbank($genome);
foreach my $i (0 .. $#lines) {
    local $_ = $lines[$i];
    next unless /^\s+CDS/;
    my ($gene_line, $gene, $desc) = infer $i;

    # CDS locations span more than one line at times, and it's up to
    # me to join them back together. Avert your eyes.
    my $expr = join ' ', @lines[$i .. $gene_line - 1];
    $expr =~ s/\r?\n\s+//g;
    my ($locs) = $expr =~ /CDS\s+(.*?)$/;
    my $seq = $parser->parse($locs);
    print "Unable to obtain $gene$/" unless $seq;

    # I extract the integers into an array by removing everything that
    # isn't. Splitting creates empty strings, which I remove.
    $locs =~ s/[^\d]/:/g;
    my @locs = grep {$_} split /:+/, $locs;

    my $start = $parser->{complemented} ? max @locs : min @locs - 12;
    my $leader = substr($genome, $start, 12);
    $leader = $parser->complement($leader) if $parser->{complemented};

    open(my $seq_h, '>' . File::Spec->catfile($dir, "$gene.txt"));
    print $seq_h "> $desc$/";
    print $seq_h $leader, $/, $seq, $/;
}
