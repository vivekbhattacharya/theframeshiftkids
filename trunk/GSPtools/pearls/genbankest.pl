use LWP::Simple qw(get);
use File::Path; use 5.010;
use strict; use warnings;

my $genomelist = $ARGV[0];

opendir(DIR, $genomelist);
my @files = grep(/\.txt$/,readdir(DIR));
closedir(DIR);

foreach my $file (@files) {
    my $accn = substr($file, 0, length($file) - 4);
    my $url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&val=%s&sendto=t";
    my $dir = $genomelist . "/" . $accn;
    
    if (!-d $dir) {
        mkpath($dir);
    }
    
    chdir($genomelist);
    
    my $s = get sprintf($url, $accn);
    open(my $h, ">$accn" . "-genbank.txt") or die "Could not download $accn.";
    print $h $s;
    close($h);
    
    system('perl C:\up\tools\pearls\genbanker.pl' . " $accn" . "-genbank.txt " . $file . ' ' . $dir);

  
}
