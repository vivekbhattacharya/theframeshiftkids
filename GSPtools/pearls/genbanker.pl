package Genbanker;
use strict; use warnings;
use File::Basename;
use lib dirname(__FILE__);
use 5.010;

use ParseGenbank;
use Smooth;

use File::Path qw(mkpath);
use File::Spec;
use List::Util qw(max min);

sub infer {
    my ($i, $big_lines) = @_;

    # Because of so many errors, we're going to place well-defined
    # bounds on how `infer` searches.
    my $max = $i + 1;
    $max++ until $max == $#{$big_lines} or $big_lines->[$max] =~ /^\s+CDS/;
    my @lines = @{$big_lines}[$i .. $max];

    my ($gene, $gene_line, $desc);
    for (my $j = 0; $j < @lines; $j++) {
        if ($lines[$j] =~ /(?:gene|locus_tag)="(.*?)"/) {
            ($gene, $gene_line) = ($1, $j);
        }
        elsif ($lines[$j] =~ /(product=|pseudo)/) {
            # Genes can have a /product= xor a /pseudo line.
            if ($1 eq 'pseudo') {$desc = '/pseudo';}
            else {
                for (join '::', @lines) {
                    s/::\s+//g; s/(\r|\n)/ /g;
                    ($desc) = /product="(.*?)"/;
                }
            }
        }
    }
    die "No gene or locus tag found starting from line $i" unless $gene;
    $desc = 'No description found' unless $desc;

    return ($gene_line + $i, $gene, $desc);
}

sub genome {
    my ($genome_file) = @_;
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
    open(my $h, ">$genome_file")
      or die 'Unable to open genome for writing.';
    print $h $genome;
    close($h);

    return $genome;
}

sub locs {
    # CDS locations span more than one line at times, and it's up to
    # me to join them back together. Avert your eyes. Precondition: $i
    # marks the beginning of a CDS line.
    my ($i, $gene_line, $lines) = @_;
    my $expr = join ' ', @{$lines}[$i .. $gene_line - 1];
    $expr =~ s/\s+//g;
    my @tokens = split m{/}, $expr;

    $expr = $tokens[0];
    $expr =~ /CDS(.*?)$/;

    unless ($1) {
        die "No locus found within the lines [$i, $gene_line]";
    }
    return $1;
}

sub parse {
    my ($genbank_file, $genome_file) = @_;
    my $genome = genome $genome_file;
    my $parser = new ParseGenbank($genome);
    my @lines = do {
        open(my $h, $genbank_file)
          or die "Could not open $genbank_file";
        <$h>;
    };

    my $i = 0;
    return sub {
        until ($lines[$i] =~ /^\s+CDS/ and $lines[$i] =~ /\d/) {
            $i += 1;
            return if $i == $#lines;
        }

        my ($gene_line, $gene, $desc) = infer($i, \@lines);
        my $locs = locs($i, $gene_line, \@lines);
        my $seq = $parser->parse($locs) or
          say "Unable to obtain $gene";

        # I extract the integers into an array by removing everything that
        # isn't. Splitting creates empty strings, which I remove.
        my @locs = grep {$_} split /[^\d]+/, $locs;
        my $start = $parser->{complemented} ? max @locs : min(@locs) - 12;
        my $leader = substr($genome, $start, 12);
        if ($parser->{complemented}) {
            $leader = Smooth::reverse_complement($leader);
        }

        $i++;
        return ($gene, $desc, $leader, $seq);
    }
}

if ($0 eq __FILE__) {
    my $genbank_file = shift or die 'No Genbank loci file given.';
    my $genome_file = shift or die 'No Genbank genome file given';
    my $dir = shift or die 'No output directory given.';

    unless (-d $dir) {
        mkpath $dir or die 'Unable to create output directory.';
    }

    my $it = parse $genbank_file, $genome_file;
    while (my ($gene, $desc, $leader, $seq) = $it->()) {
        open(my $h, '>', File::Spec->catfile($dir, "$gene.txt"));
        say $h "> $desc";
        say $h $leader;
        say $h $seq;
    }
}
1;
