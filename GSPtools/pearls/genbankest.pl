use LWP::Simple qw(get);
use File::Path; use 5.010;
use File::Spec;
use strict; use warnings;

chdir($ARGV[0]);
my @files = glob('*.txt');

my $url = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&val=%s&sendto=t';
foreach my $file (@files) {
    my $accn = substr($file, 0, -4);
    my $dir = File::Spec->catdir($genomelist, $accn);

    mkpath $dir if !-d $dir;
    chdir($genomelist);

    my $s = get sprintf($url, $accn);
    open(my $h, ">$accn-genbank.txt") or die "Could not download $accn";
    print $h $s;

    system(sprintf('perl "C:\up\tools\pearls\genbanker.pl %s-genbank.txt %s %s', $accn, $file, $dir));
}
