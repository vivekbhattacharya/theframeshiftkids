#!/usr/bin/env perl
use strict; use warnings; use 5.010;
use File::Basename;
use Test::Simple tests => 2;

sub check {
    my ($seq, $p, $standard) = @_;

    my $candy = `perl ../scan_brightly.pl -p $p auuccuccacuag $seq`;

    open(my $handle, '<', $standard);
    my $cavities = <$handle>;
    $cavities =~ s/\s*//g;
    $candy =~ s/\s*//g;

    ok($candy eq $cavities, "$seq with $p");
}

check 'prayforme.txt', 'Kidnap::Freier', 'prfB-Freier.txt';
check 'prayforme.txt', 'Kidnap::XiaMathews', 'prfB-XiaMathews.txt';
