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
    while (1) {
        last if $max >= $#{$big_lines};
        last if $big_lines->[$max] =~ /^\s+(CDS|gene)/;
        $max++;
    }
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

        # Thank you, NC_007456's DNA polymerase.
        unless ($gene) {
            say "No name found (near line $i, makes $desc), defaulting to \"unknown\"";
            $gene = "unknown";
        }

        # Thank you, NC_007817 and your A* gene.
        $gene =~ s/\*/-star/g;

        # Thank you, NC_007204 and your lysS/U genes.
        $gene =~ s{/}{-}g;
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
        if ($seq =~ /[^aucg]/) {
            say "$gene has a non-AUCG base.";
        }

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

if ($0 eq __FILE__) {
    Smooth::helpcheck();
    main @ARGV;
}

1;

__END__

=head1 NAME

genbanker.pl

=head1 SYNOPSIS

genbanker.pl I<Genbank Nucleotide> I<Genbank Genome> I<output
directory>

=head1 DESCRIPTION

Genbanker extracts all the genes from a Genbank Genome entry given the
corresponding Genbank Nucleotide entry. To grab a Genome, I recommend
F<pull_genome.pl>. To grab a Nucleotide and automate the call to
Genbanker, I recommend F<genbankeset.pl>.

=head1 OPTIONS AND ARGUMENTS

If passed with no arguments, this help text pops up.

=over

=item I<Genbank Nucleotide>

The easiest way to find the Nucleotide file is automatically via
F<genbankest.pl>. However, you can also search the Nucleotide database
online. The Genbank website's search engine is subpar, so unless you
know the accession number for your genome, this might be difficult.
Once you find it, format it as plain text by clicking "Send to" and
selecting "Text". Save this somewhere. Pass the location here.

[1]: http://www.ncbi.nlm.nih.gov/entrez?db=nucleotide

=item I<Genbank Genome>

The easiest way to find the Genome file is automatically via
F<pull_genome.pl>. You can also search Genome database online [1]. For
an index of genomes, consult the sidebar. For example, the Bacteria
index [2]. Again, send this to a text file and pass the path to me.

=over

=item [1]: http://www.ncbi.nlm.nih.gov/sites/entrez?db=genome

=item [2]: http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=2&type=1&name=Bacteria%20Complete%20Chromosomes

=back

=item I<output directory>

This is where I will save all the genes. I will the create the
directory if it does not already exist. I save each gene to a FASTA
file named "gene name.txt" where gene name is the one listed in the
CDS subsection for that gene in the Nucleotide file. The listing
occurs in /gene= or /locus_tag=. If neither exists, I spit out a
message and call it "unknown" instead. The FASTA description is pulled
from /product= or /pseudo or /note= in that order. Otherwise, it's "No
description found."

I automatically pull leader sequences and complement mRNA sequences. I
know to pull leader sequences after complemented genes. I'm just that
good. Gene names with asterisks in them, like A*, are dealt by
renaming them to A-star. Duplicate gene names are resolved by
appending "-again" and then a random integer.

I rarely fail, but when I do make a mistake I do it fabulously. Dance
like nobody is watching.

=back

=head1 COPYRIGHT AND LICENSE

I'm under BSD. See the License article on the Google Code wiki.
