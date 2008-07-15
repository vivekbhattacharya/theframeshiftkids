package Genbankest;
use File::Basename;
use lib dirname(__FILE__);

package Genbanker;
require 'genbanker.pl';

package Genbankest;
use LWP::Simple qw(get);
use File::Path qw(mkpath);
use File::Spec;
use strict; use warnings; use 5.010;

chdir($ARGV[0]);
my @files = glob('*.txt');

my $url = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&list_uids=%s&dopt=gbwithparts&sendto=t&fmt_mask=295416';
foreach my $file (@files) {
    next if $file =~ /-genbankest\.txt/;

    my $accn = substr($file, 0, -4);
    my $dir = File::Spec->catdir($ARGV[0], $accn);
    mkpath $dir if -d $dir;

    my $s = get sprintf($url, $accn);
    die "Empty Nucleotide web page for $file" unless $s =~ /\S/;

    open(my $h, '>', "$accn-genbankest.txt") or die "Could not download $accn";
    print $h $s;

    say "$accn-genbankest.txt, $dir, $file";
    Genbanker::main("$accn-genbankest.txt", $file, $dir);
}
