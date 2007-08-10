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

use File::Path qw(mkpath);
use Smooth;
if (__FILE__ eq $0) {
    Smooth::helpcheck();
    my ($file, $folder) = @ARGV;
    open(my $genehandle, $file) or die('Genes.');
    
    mkpath($folder);
    while (my $gene = <$genehandle>) {
        chomp $gene;
        my $id = get_id $gene;
        
        if ($id =~ m|[\s/]|) {
            print "Could not find gene $gene", $/;
            next;
        }
        print $gene, $/;
        
        open(my $handle, ">$folder/$gene.txt");
        print $handle thehardpart $id;
    }
}

1;

__END__

=head1 NAME

happy_face.pl

=head1 SYNOPSIS

    happy_face.pl genes.txt c:\output

=over 20

=item B<happy_face.pl>

I<list of genes> I<output path>

=back

=head1 DESCRIPTION

happy_face screenscrapes the EcoGene database website for
genes listed in a file, separated by newlines. It then saves
the gene sequence to "gene.txt" in the output path, where
"gene" is actually the gene sequence name, such as aceF or
boobies.

=head1 CAVEATS

Liable to broken as soon as EcoGene updates its website.