use Smooth;
use warnings;
use strict;

sub printo {
    my $codons = shift;
    my @codons = @$codons;
    
    my $i = 0;
    my $total = scalar @codons;
    while(1) {
        my $end = $i + 9;
        # In the common case the last codon doesn't
        # end on a multiple of 10, Perl will slice
        # some undefs for us.
        $end = $total % 10 + $i - 1 if $end >= $total;
        
        my @subcodons = @codons[$i .. $end];
        print join ' ', @subcodons, $/;
        
        $i += 10;
        last if $i > @codons;
    }
}

sub run {
    my @codons = ();
    Smooth::webopen shift @ARGV, sub {
        my $line = $_;
        
        # http://en.wikipedia.org/wiki/List_of_standard_amino_acids#Gene_expression_and_biochemistry
        # Sanitize and check for empty lines.
        $line =~ s/[^A-Z]//g;
        my @chars = split //, $line;
        return unless $line;
        
        map { push @codons, Smooth::cupid($_); } @chars;
    };
    push @codons, 'uga';
    print 'gcc aua ggc uau', $/;
    printo(\@codons);
}

sub check {
    ## Get a list of codons.
    my @codons = ();
    Smooth::webopen shift @ARGV, sub {
        my $line = $_;
        push @codons, split(/\s+/, $line);
    };
    # -2 because the last codon is the stop codon.
    @codons = @codons[4 .. @codons-2];
    
    ## Convert that list to proteins.
    my $proteins = '';
    map {
        $proteins .= Smooth::reverse_cupid($_);
    } @codons;
    
    ## Get actual list of proteins.
    my $actual = '';
    Smooth::webget shift @ARGV, sub {
        $actual = join '', @_;
    };
    $actual =~ s/[^A-Z]//g;
    
    if ($actual eq $proteins) { print "# They are equal. Congratulations.", $/; }
    else { print "# They are not equal.", $/; }
    
    print "Codons: $actual$/";
    print "Actual: $proteins$/";
}

if ($0 eq __FILE__) {
    my $help = <<END;
NAME
    dysentery.pl

USAGE
    1) Convert protein sequence (1-letter abbreviation standard) to codons.
        `perl dysentery.pl [url/path to protein sequence]`
    2) Check if a codon sequence produces a protein sequence.
        `perl dysentery.pl [url/path to codons] [url/path to proteins]`
END
    if (!@ARGV or $ARGV[0] eq '--help') { print $help; }
    elsif ($ARGV[0] eq '--check') { shift @ARGV; check(); }
    else { run(); }
}
