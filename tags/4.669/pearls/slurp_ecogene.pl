use strict; use warnings;
use Text::Wrap;
use Smooth;
$Text::Wrap::columns = 60;

Smooth::helpcheck();

my ($genome_file, $points_file) = @ARGV;
my $genome = lc Smooth::webslurp $genome_file;
$genome =~ s/\n//g;
$genome =~ tr/t/u/;

sub extraction {
    my ($start, $end) = @_;
    return substr($genome, $start - 1, $end - $start + 1);
}

mkdir('genes'); mkdir('confused');
Smooth::webopen $points_file, sub {
    chomp;
    my ($gene, $start, $end, $desc) = split "\t";
    my $complement = 0;
    if ($start < 0) {
        $start = -$start;
        $complement = 1;
    }

    my $seq = extraction($start, $end);
    my $leader = '';
    my $confused = 0;
    if ($seq =~ /^(aug|gug)/) {
        # Oh good, 5' to 3'.
        $leader = extraction($start - 12, $start - 1);
    }
    else {
        if ($complement or $seq =~ /(cau|cac)$/) {
            # What the hell? 3' to 5'?
            $leader = extraction($end + 1, $end + 12);
            $seq = scalar reverse $seq;
            $seq =~ tr/aucg/uagc/;

            $leader = scalar reverse $leader;
            $leader =~ tr/aucg/uagc/;
        }
        else {
            $confused = 1;
        }
    }
    $seq = wrap('', '', $seq);
    if ($confused) {
        open(my $g, ">confused/$gene.txt");
        print $g ">$desc (slurp: orientation unknown)\n$seq";
    }
    else {
        open(my $f, ">genes/$gene.txt");
        print $f ">$desc\n$leader\n$seq";
    }
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

It automatically detects whether the gene is written 3' to 5' in the
genome or otherwise, albeity crudely. Genes that do not fit this status
quo are placed in a "confused" folder.

=head1 SYNOPSIS

=over 20

=item B<slurp_ecogene.pl>

[B<--help>]

=item B<slurp_ecogene.pl>

I<genome file> I<points file>

=back

=cut
