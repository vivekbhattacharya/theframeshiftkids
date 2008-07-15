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
    $max++ until $max == $#{$big_lines} or $big_lines->[$max] =~ /^\s+(CDS|gene)/;
    my @lines = @{$big_lines}[$i .. $max];

    for (join ' ', @lines) {
        s/(\r|\n)/::/g; s/(\s+|::)/ /g;
        my ($locs) = /CDS\s+(.*?)\//;
        my ($gene) = /(?:gene|locus_tag)="(.*?)"/;
        my ($desc) = /product="(.*?)"/;

        if (!$desc) {
            $desc = 'No description found';
            if (/pseudo/) {$desc = '/pseudo';}
            if (/note="(.*?)"/) {$desc = "/note = $1";}
        }

        # Not sure why there's a '<' in front of some genes.
        # Addendum: Sometimes, there's a '>'! (NC_005841)
        $locs =~ s/[<>]//g;

        # Thank you NC_007456's DNA polymerase.
        unless ($gene) {
            say "No name found (near line $i, makes $desc), defaulting to \"unknown\"";
            $gene = "unknown";
        }

        # Thank you NC_007817 and your A* gene.
        $gene =~ s/\*/-star/g;
        return ($gene, $desc, $locs);
    }
}

sub genome {
    my ($genome_file) = @_;
    my $genome = lc do {
        local $/;
        open(my $genome_h, '<', $genome_file)
          or die "Unable to open genome $genome_file.";
        <$genome_h>;
    };

    $genome =~ s/^>.*$//m; # Get rid of FASTA header.
    $genome =~ s/\s//g; # Get rid of newlines.
    $genome =~ tr/t/u/; # Convert to mRNA.

    # Write back a sanitized version for faster future runs.
    open(my $h, '>', $genome_file)
      or die 'Unable to open genome for writing.';
    print $h $genome;
    close($h);

    return $genome;
}

sub parse {
    my ($genbank_file, $genome_file) = @_;
    my $genome = genome $genome_file;
    my $parser = new ParseGenbank($genome);
    my @lines = do {
        open(my $h, $genbank_file)
          or die "Could not open Genbank Nucleotide file $genbank_file";
        <$h>;
    };

    my $i = 0;
    return sub {
        # Bad hack: Some amino acid lines have 'CDS'. Testing for a
        # digit should produce the correct result.
        until ($lines[$i] =~ /^\s+CDS.*?\d/) {
            $i += 1;
            return if $i == $#lines;
        }

        my ($gene, $desc, $locs) = infer($i, \@lines);
        my $seq = $parser->parse($locs) or say "Unable to obtain $gene";

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

sub main {
    my $genbank_file = shift or die 'No Genbank loci file given.';
    my $genome_file = shift or die 'No Genbank genome file given';
    my $dir = shift or die 'No output directory given.';

    unless (-d $dir) {
        mkpath $dir or die 'Unable to create output directory.';
    }

    my $it = parse $genbank_file, $genome_file;

    my %already = ();
    while (my ($gene, $desc, $leader, $seq) = $it->()) {
        # Loop until we find a unique filename. -f won't work because
        # it's possible the output directory is dirty from the last
        # genbanker.pl run.
        while (1) {
            if (!$already{$gene}) {
                $already{$gene} = 1; last;
            }
            my $x = int(rand(100));
            $gene .= "-again-$x";
        }
        my $file = File::Spec->catfile($dir, "$gene.txt");
        open(my $h, '>', $file) or die "Unable to write to $file";
        say $h "> $desc\n$leader\n$seq";
    }
}

main @ARGV if $0 eq __FILE__;

1;
