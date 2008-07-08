# All the functions we'll never call from the command
# line but need anyway.
package Smooth;
use LWP::Simple qw(get);
use strict; use warnings;
BEGIN {
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT_OK = qw(prot2codon codon2prot getseq);
}

# Maps $func to the lines of the $file, be it on the Internet or not.
# Obscure trivia: Like map, it returns the resulting array but nobody
# uses it except getseq.
sub webopen {
    my ($file, $func) = @_;
    if ($file =~ m|^http://|) {
        my @lines = split($/, get $file);
        map { $func->($_) } @lines;
    } else {
        open(my $handle, $file) or die "webopen: Cannot open file `$file`";
        map { $func->($_) } <$handle>;
    }
}

# Read the file/URL and return the contents.
sub webslurp {
    my ($file) = @_;
    if ($file =~ m|^http://|) {
        return get $file;
    } else {
        open(my $handle, $file) or die "webslurp: Cannot open file `$file`";
        local $/; return <$file>;
    }
}

sub getseq {
    # FASTA files are so annoying.
    my @lines = webopen shift, sub { $_ unless /^>/ };
    sanitize(join '', @lines);
}

sub seq2codons {
    my $seq = shift;
    # Split the sequence into threes using a nasty
    # trick involving exclamation points and intrigue.
    $seq =~ s/(...)/$1!/g;
    split '!', $seq;
}

# Convert all gene sequences to this normalized form.
sub sanitize {
    local $_ = shift;
    s/[\s0-9]//g; $_ = lc $_; tr/t/u/; chomp; $_;
}

# Cf. page 43-6 of _Programming Perl for Bioinformatics_
# by James D. Tisdall
sub reverse_complement {
    for(shift) { tr/autcg/uaagc/; return scalar reverse; }
}

# http://en.wikipedia.org/wiki/List_of_standard_amino_acids
our %expression = (
    'A' => 'gcu,gcc,gca,gcg',
    'C' => 'ugu,ugc',
    'D' => 'gau,gac',
    'E' => 'gaa,gag',
    'F' => 'uuu,uuc',
    'G' => 'ggu,ggg,gga,ggc',
    'H' => 'cau,cac',
    'I' => 'auu,auc,aua',
    'K' => 'aaa,aag',
    'L' => 'uua,uug,cuu,cua,cuc,cug',
    'M' => 'aug',
    'N' => 'aau,aac',
    'P' => 'ccu,cca,ccg,ccc',
    'Q' => 'caa,cag',
    'R' => 'cga,cgu,cgc,cgg,aga,agg',
    'S' => 'ucu,uca,ucg,ucc,agu,agc',
    'T' => 'acu,acc,aca,acg',
    'V' => 'guu,guc,gua,gug',
    'W' => 'ugg',
    'Y' => 'uau,uac',
    # Stop codons
    '.' => 'uga,uaa,uag',
);

# Generate a reverse hashmap forthwith!
our %repression = ();
while (my ($key, $value) = each(%expression)) {
    my @codons = split /,/, $value;
    map { $repression{$_} = $key; } @codons;
}

# Randomly pulls the equivalent codon given a
# one-character uppercase
#
# Test protein sequence: http://shadytrees.pastebin.ca/raw/551240
# Test gene sequence: http://shadytrees.pastebin.ca/raw/551246
sub prot2codon {
    map {
        my @codons = split /,/, $expression{$_};
        # Scalar context goodness
        $codons[rand @codons];
    } @_;
}

sub codon2prot {
    map { $repression{$_} } @_;
}

use Pod::Usage;
sub helpcheck {
    if (!@ARGV or $ARGV[0] eq '--help') {
        pod2usage(-verbose => 3, -noperldoc => 1); exit;
    }
}
1;
