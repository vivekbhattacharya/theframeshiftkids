package Hansel;
use strict; use warnings;
use lib '../pearls';

use Smooth qw(getseq);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

# World's simplest diff function
sub diff {
    my ($old, $new, $start) = @_;
    my @old = @$old;
    my @new = @$new;
    
    die 'Two sequences are not of equal length' unless @old == @new;
    
    for my $i ($start + 1 .. @old) {
        next unless $old[$i] && $new[$i];
        next if $old[$i] eq $new[$i];
        return $i;
    }
    return -1;
}

sub backtrack {
    my ($old, $new, $start, $folder) = @_;
    my @old = Smooth::seq2codons getseq $old;
    my @new = Smooth::seq2codons getseq $new;
    my $i = diff \@old, \@new, $start;

    my $janet = join $/, @old;
    
    $old[$i] = $new[$i];
    my $reno = join $/, @old;
    
    ($janet, $reno, $i, $folder);
}

use File::Path;
sub rite {
    my ($old, $new, $i, $folder) = @_;
    $folder .= "/$i";
    mkpath $folder;
    
    if ($i == -1) {
        print "{-1 -1}";
        exit;
    }
    
    open my $cat, ">$folder/old.txt";
    print $cat $old;
    
    open my $dead_cat_store,  ">$folder/new.txt";
    print $dead_cat_store $new;
    
    # Talk with Matlab
    print "{$i '$folder'}";
}

sub print_diff {
    my ($old, $new, $start, $folder) = @_;
    my @old = Smooth::seq2codons getseq $old;
    my @new = Smooth::seq2codons getseq $new;
    
    print "# Index, old codon, new codon", $/;
    for (
        my $i = diff(\@old, \@new, -1); $i > -1;
        $i = diff(\@old, \@new, $i)
    ) {
        print "$i, $old[$i], $new[$i]", $/;
    }
    
}

if ($0 eq __FILE__) {
    my ($help, $diff, $backtrack);
    GetOptions('help|?' => \$help, 'diff' => \$diff,
        'backtrack' => \$backtrack) or pod2usage(-verbose => 3);
    pod2usage(-verbose => 3) if ($help || !@ARGV);
    
    if ($diff) { print_diff @ARGV }
    elsif ($backtrack) { rite backtrack @ARGV }
    else { pod2usage(-verbose => 1) }
}

1;

__END__

=head1 NAME

hansel.pl

=head1 SYNOPSIS

    hansel.pl --backtrack old-rpoS.txt new-rpoS.txt 46 "C:\Work folder"

=cut