package HappyFace;
use strict;
use warnings;
use LWP::Simple;

my ($file) = @ARGV;

open(my $genehandle, $file) or die('Genes.');

my $first = 'window.location = ("geneinfo.php?eg_id=';

while (my $gene = <$genehandle>) {
    chomp $gene;
    my $wholeurl = "http://ecogene.org/ecoSearchProcess.php?name=$gene&numofresults=20&searchType=gene";
    my $site = get($wholeurl);
    my $begin = index $site, $first;
    $begin = $begin + length($first);
    my $end = index $site, '"', $begin;
    my $eg = substr $site, $begin, $end - $begin;
    print $eg;
}

1;