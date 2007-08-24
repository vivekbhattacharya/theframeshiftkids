use strict; use warnings;
use File::Basename;
use Test::Simple tests => 2;

sub check {
    my ($seq, $p, $standard) = @_;
    
    my $candy = `perl ..\\scan_brightly.pl -p $p auuccuccacuag $seq`;
    local $/;
    open(my $handle, $standard);
    my $cavities = <$handle>;
    ok($candy eq $cavities, "$seq with $p");
}

check 'prfB.txt', 'Kidnap::Freier', 'prfB-Freier.txt';
check 'prfB.txt', 'Kidnap::XiaMathews', 'prfB-XiaMathews.txt';