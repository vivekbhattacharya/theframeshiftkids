# All the functions we'll never call from the command
# line but need anyway.
package Smooth;
use LWP::Simple qw(get);

sub webopen {
    my ($file, $func) = @_;
    if ($file =~ m|^http://|) {
        map &$func, split($/, get($file));
    } else {
        open(my $handle, $file) or die "webopen: Cannot open file `$file`";
        map &$func, <$handle>;
    }
}

1;