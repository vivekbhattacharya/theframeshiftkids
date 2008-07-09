use LWP::Simple qw(get);
use File::Path; use 5.010;
use File::Spec;
use strict; use warnings;

chdir($ARGV[0]);
my @files = glob('*.txt');

my $url = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&val=%s&sendto=t';
foreach my $file (@files) {
    next if $file =~ /-genbankest\.txt/;

    my $accn = substr($file, 0, -4);
    my $dir = File::Spec->catdir($ARGV[0], $accn);
    mkpath $dir if -d $dir;

    my $s = get sprintf($url, $accn);
    die "Empty genome $file" unless $s =~ /\S/;

    open(my $h, '>', "$accn-genbankest.txt") or die "Could not download $accn";
    print $h $s;

    my $cmd = sprintf('perl "C:\up\tools\pearls\genbanker.pl" %s-genbankest.txt %s "%s"', $accn, $file, $dir);
    say $cmd;
    system($cmd);
}
