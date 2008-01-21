use strict; use warnings;

unless (@_) {
    print 'I need the absolute path to genbank2ecogene.pl.', $/;
    exit;
}

my $temp = 'thisisrubbishiapologize';
my $script = shift;

`perl $script goody.txt > $temp`;
my $output = `diff -u points.txt $temp`;

if ($output) {
    print 'FAIL', $/, $output, $/;
}
else {
    print 'ok', $/;
}

unlink($temp);
