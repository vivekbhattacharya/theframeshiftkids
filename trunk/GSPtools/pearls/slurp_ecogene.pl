use strict;
use warnings;
use Text::Wrap;
use Smooth;
$Text::Wrap::columns = 60;

Smooth::helpcheck();

my ($genome_file, $points_file) = @ARGV;
my $genome = lc Smooth::webslurp $genome_file;
$genome =~ s/\n//g;
$genome =~ tr/atgc/uacg/;

mkdir('genes');
my %points;
Smooth::webopen $points_file, sub {
    my ($gene, $start, $end, $desc) = split "\t";
    open(my $f, ">genes/$gene.txt");

    my $leader = substr($genome, $end, 12);
    my $seq = substr($genome, $start - 1, $end - $start + 1);

    $leader = scalar reverse $leader;
    $seq = wrap('', '', scalar reverse $seq);
    print $f ">$desc\n$leader\n$seq";
};
1;

__END__

=head1 NAME

ecogene_slurp.pl (web-enabled): Reads a genome and Ecogene
tab-delimited points file to extract gene sequences into the "genes"
folder. See L<http://ecogene.org/DatabaseTable.php> for a points file.
See L<http://ecogene.org/SequenceDownload.php> for the genome. You
must select the headers "Gene Name", "Left End", "Right End", and
"Description" and no others.

=head1 SYNOPSIS

=over 20

=item B<cornerstone.pl>

[B<--help>]

=item B<cornerstone.pl>

I<genome file> I<points file>

=back

=cut
