package HappyFace;
use strict;
use warnings;
use LWP::Simple;

my ($file) = @ARGV;

open(my $genehandle, $file) or die('Genes.');

my $begin_ecogene = 'http://ecogene.org/ecoSearchProcess.php?name=';
my $end_ecogene = '&numofresults=20&searchType=gene';
my $first = 'window.location = ("geneinfo.php?eg_id=';

while (my $gene = <$genehandle>) {
    chomp $gene;
    my $wholeurl = $begin_ecogene . $gene . $end_ecogene;
    my $site = get($wholeurl);
    my $begin = index $site, $first;
    $begin = $begin + length($first);
    my $end = index $site, '"', $begin;
    my $eg = substr $site, $begin, $end - $begin;
    print $eg;
}

1;