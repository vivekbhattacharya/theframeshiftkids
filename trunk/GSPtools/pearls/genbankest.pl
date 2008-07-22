package Genbankest;
use File::Basename;
use lib dirname(__FILE__);

package Genbanker;
require 'genbanker.pl';

package Genbankest;
use Getopt::Long;
use LWP::Simple qw(get);
use File::Path qw(mkpath);
use File::Spec;
use strict; use warnings; use 5.010;

my $download = 1;
my $help = 0;
GetOptions('download!' => \$download, 'help|?' => \$help) or Smooth::help();
Smooth::help() if $help or not @ARGV;

chdir($ARGV[0]);

my $url = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&list_uids=%s&dopt=gbwithparts&sendto=t&fmt_mask=295416';
foreach my $file (glob('*.txt')) {
    next if $file =~ /-genbankest\.txt/;

    # Strip the .txt off.
    my $accn = substr($file, 0, -4);
    my $dir = File::Spec->catdir('genbankest', $accn);
    mkpath $dir if -d $dir;

    my $nucleotide = "$accn-genbankest.txt";
    unless (-f $nucleotide and $download) {
        my $s = get sprintf($url, $accn);
        die "Empty Nucleotide web page for $file" unless $s =~ /\S/;
        open(my $h, '>', $nucleotide) or die "Could not download $accn";
        print $h $s;
    }

    say "$accn-genbankest.txt, $dir, $file";
    Genbanker::main("$accn-genbankest.txt", $file, $dir);
}

__END__

=head1 NAME

genbankest.pl

=head1 SYNOPSIS

genbankest.pl [B<--help>] [B<--(no)-download>] I<directory of genomes>

=head1 DESCRIPTION

For each genome in the directory, Genbankest downloads the necessary
Nucleotide entry from Genbank and runs the pair against Genbanker.
This quickly creates directories of all the possible pulled genes from
all the genomes.

=head1 OPTIONS AND ARGUMENTS

All options can be abbreviated thanks to Getopt::Long.

=over

=item I<directory of genomes>

Each genome should be of the format "[accession number].txt". Each
file should contain an mRNA sequence in a plain-text or FASTA
plain-text format.

=item --help

Increases rainfall by an improbable likelihood.

=item --no-download, --download

The default, if unspecified, is --no-download. Genbankest will not
download a Nucleotide entry if one already exists. Use --download to
flush this rudimentary cache.

=back

=head1 CAVEATS

The Genbank file formats are unpredictable, messy, and poorly
designed. Genbanker and Genbankest do their best to catch all the
corner cases, but sometimes they explode quite spectacturly (with
lights and tinsel and I<everything>) when they cannot. You know who to
email if that happens.

=head1 COPYRIGHT AND LICENSE

I'm under BSD. See the License article on the Google Code wiki.
