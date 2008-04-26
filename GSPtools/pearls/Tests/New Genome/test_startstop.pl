use strict;
use warnings;
use File::Basename;
use lib dirname(__FILE__) . '/../..';
use Smooth qw(getseq);

chdir(shift @ARGV);
foreach (glob('*.txt')) {
    my $a = Smooth::getseq $_;
    $a = substr($a, 12);
    if ($a =~ /^(gug|aug)/) {next;}
    else {print "$_ does not have a start codon\n";}

    if ($a =~ /(uag|uga|uaa)$/) {next;}
    else {print "$_ does not have a stop codon\n";}
}
