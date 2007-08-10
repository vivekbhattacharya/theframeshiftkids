package HappyFace;
use strict; use warnings;
use LWP::Simple;

our $tryst = 'window.location = ("geneinfo.php?eg_id=';
sub get_id {
    my ($gene) = @_;
    my $wholeurl = "http://ecogene.org/ecoSearchProcess.php?name=$gene&numofresults=20&searchType=gene";
    my $site = get($wholeurl);

    my $begin = index($site, $tryst) + length($tryst);
    my $end = index $site, '"', $begin;
    return substr $site, $begin, $end - $begin;
}

sub thehardpart {
    my ($id) = @_;
    my $url = "http://ecogene.org/dnaSequence.php?eg_id=$id&jointext=&us=12&ds=0";
    my $site = get($url);
    
    my $begin = index $site, ">$id";
    my $end = index $site, '</PRE>', $begin;
    return substr $site, $begin, $end - $begin;
}


if (__FILE__ eq $0) {
    my ($file) = @ARGV;
    open(my $genehandle, $file) or die('Genes.');

    while (my $gene = <$genehandle>) {
        chomp $gene;
        my $id = get_id $gene;
        print thehardpart $id;
    }
}

1;