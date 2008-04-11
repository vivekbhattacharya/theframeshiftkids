use strict;
use warnings;
use lib '../..';
use Smooth qw(getseq);

chdir('genes');
foreach (glob('*.txt')) {
    my $a = Smooth::getseq $_;
    $a = substr($a, 12);
    if ($a =~ /^(gug|aug)/) {next;}
    else {print '$_ does not have a start codon';}

    if ($a =~ /(uag|uga|uaa)$/) {next;}
    else {print '$_ does not have a stop codon';}
    print "\n";
}
