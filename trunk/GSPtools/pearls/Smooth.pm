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

# Yields an array of lines from $file to $func,
# be it on the Internet or not. The array
# is NOT a reference.
sub mooch {
    my ($file, $func) = @_;
    my @lines = ();
    
    if ($file =~ m|^http://|) {
        my $contents = get $file;
        @lines = split($/, $contents);
    } else {
        open(my $handle, $file) or die "webget: Cannot open file `$file`";
        @lines = <$handle>;
    }
    # Handle Windows line endings.
    map { s/\r|\n//g } @lines;
    &$func(@lines);
}

# Maps $func to the lines of the $file,
# be it on the Internet or not.
sub webopen {
    my ($file, $func) = @_;
    # The following rates 9/10 on the Awesome Scale.
    mooch $file, sub {
        map { $func->($_) } @_;
    };
}

# Read the file/url and return the contents.
sub webslurp {
    mooch shift, sub {
        return join('', @_);
    };
}

sub getseq {
    my $file = shift;
    mooch $file, sub {
        for(join '', @_) { return sanitize($_) }
    };
}

# Convert all gene sequences to this normalized form.
sub sanitize {
    $_ = shift;
    s/[\s0-9]//g; tr/A-Z/a-z/; tr/t/u/; $_;
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

sub permute {
    use List::Util qw(reduce);
    no warnings 'once';
    # Note: we have to store a temporary
    # $x to access $_ after it's been shadowed.
    reduce {[
        map {
            my $x = $_;
            map "$_  $x", @$a;
        } @$b
    ]} @_;
}
1;