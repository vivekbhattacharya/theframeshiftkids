use strict; use warnings;
use lib '../pearls';

package Hansel;
use Smooth qw(getseq);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Temp qw(tempfile);

sub diff {
    my ($a, $b) = @_;
    my @janet = Smooth::seq2codons getseq $a;
    my @reno = Smooth::seq2codons getseq $b;
    
    die 'Two sequences are not of equal length' unless @janet == @reno;
    
    my @changes;
    for my $i (0 .. @janet) {
        next unless $janet[$i] && $reno[$i];
        next if $janet[$i] eq $reno[$i];
        push @changes, "$i,$janet[$i],$reno[$i]";
    }
    
    my $str = join $/, @changes;
    my $handle = tempfile();
    print $handle $str;
    $str;
}

sub print_diff {
    print Dumper @_;
}

if ($0 eq __FILE__) {
    my ($help, $diff);
    GetOptions('help|?' => \$help, 'diff' => \$diff) or pod2usage(-verbose => 3);
    pod2usage(-verbose => 3) if ($help || !@ARGV);
    
    print_diff diff shift, shift if $diff;
}

1;

__END__

=head1 NAME

hansel.pl

=cut