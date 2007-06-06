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
    webget $file, sub { map &$func, @_; };
}

sub getseq {
    my $file = shift;
    webget $file, sub {
        for(join '', @_) { s/[\s0-9]//g; tr/A-Z/a-z/; tr/t/u/; return $_; }
    };
}
1;