use strict; use warnings;
use Digest::SHA1;

unless (@ARGV) {
    die 'I need the absolute path to genbanker.pl.';
}

my $temp = 'thisisrubbishiapologize';
my $script = shift @ARGV;

`perl "$script" goody.txt genome.txt $temp`;

my $sha1 = Digest::SHA1->new;
foreach (qw/a b c/) {
    open(my $h, "$temp/$_.txt");
    $sha1->addfile($h);
    close($h);
}

use Data::Dumper;

# Be careful! hexdigest destroys the state of $sha1!
if ($sha1->hexdigest eq 'c205849564a51bfc54b8da7d793f46f411146b49') {
    print 'ok 1', $/;
}
else {
    print 'FAIL', $/;
}

use File::Basename qw(dirname);
use File::Path qw(rmtree);
use File::Spec;
rmtree(File::Spec->catdir(dirname(__FILE__), $temp));
