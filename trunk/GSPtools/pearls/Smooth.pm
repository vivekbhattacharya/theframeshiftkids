# All the functions we'll never call from the command
# line but need anyway.
package Smooth;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(prot2codon codon2prot);
use LWP::Simple qw(get);

# Yields an array of lines from $file to $func,
# be it on the Internet or not. The array
# is NOT a reference.
sub webget {
    my ($file, $func) = @_;
    if ($file =~ m|^http://|) {
        # LWP::Simple
        my $contents = get $file;
        my @lines = split($/, $contents);
        # Handle Windows line endings.
        map { s/\r|\n//g } @lines;
        &$func(@lines);
    } else {
        open(my $handle, $file) or die "webget: Cannot open file `$file`";
        my @lines = <$handle>;
        # Handle Windows line endings.
        map { s/\r|\n//g } @lines;
        &$func(@lines);
    }
}

# Maps $func to the lines of the $file,
# be it on the Internet or not.
sub webopen {
    my ($file, $func) = @_;
    # The following rates 9/10 on the Awesome Scale.
    webget $file, sub { map { &$func(shift); } @_; };
}

sub getseq {
    my $file = shift;
    webget $file, sub {
        for(join '', @_) { s/[\s0-9]//g; tr/A-Z/a-z/; tr/t/u/; return $_; }
    };
}

# http://en.wikipedia.org/wiki/List_of_standard_amino_acids
my %expression = (
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
);

# Generate a reverse hashmap forthwith!
my %repression = ();
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
    my @codons = split /,/, $expression{+shift};
    # Scalar context goodness
    $codons[rand @codons];
}

sub codon2prot { $repression{+shift} }
1;