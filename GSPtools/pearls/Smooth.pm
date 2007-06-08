# All the functions we'll never call from the command
# line but need anyway.
package Smooth;
use LWP::Simple qw(get);

# Yields an array of lines from $file to $func,
# be it on the Internet or not.
sub webget {
    my ($file, $func) = @_;
    if ($file =~ m|^http://|) {
        &$func(split($/, get $file));
    } else {
        open(my $handle, $file) or die "webget: Cannot open file `$file`";
        &$func(<$handle>);
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

# Randomly pulls the equivalent codon given a
# one-character uppercase 
sub cupid {
    my $char = shift;
    my $codons = $expression{$char};
    my @codons = split /,/, $codons;
    
    # Scalar context goodness
    $codons[rand @codons];
}
1;