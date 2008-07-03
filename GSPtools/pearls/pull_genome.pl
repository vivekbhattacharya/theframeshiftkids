use LWP::Simple qw(get);
use File::Path; use 5.010;
use strict; use warnings;

my $url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?list_uids=%s&dopt=fasta&sendto=t\n";
my $dir = $ARGV[0];
if (!-d $dir) {
    mkpath($dir);
    say "$dir created";
}
chdir($ARGV[0] or die 'No folder given.') or die 'Non-existent directory.';

while (<stdin>) {
    chomp;
    my $s = get sprintf($url, $_);
    open(my $h, ">$_.txt") or die "Could not download $_.";
    print $h $s;
    close($h);
}

__END__

=head1 NAME

get_genome.pl: Gets genome from the internet given a Genbank accession
number

=head1 SYNOPSIS

type cheese.txt | perl get_genome.pl "c:\output directory"

=head1 DESCRIPTION

I read from STDIN a list of accession numbers and saves each one to a
text file in the given output directory. The output directory does not
need to exist; I will create it if necessary. Because I'm reading from
STDIN, you can save a list of accn.'s separated by newlines and then
pipe it into get_genome.pl with Windows' `type` or *nix's `cat`. This
is demonstrated in the SYNOPSIS.

=head1 SEE ALSO

genbanker.pl: parses genome file and a Genbank loci listing to give
you genes

=head1 COPYRIGHT

This is licensed under the Beached Whale Frameshift project.

