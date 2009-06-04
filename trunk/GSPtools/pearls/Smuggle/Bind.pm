#!/usr/bin/env perl

# Based on code by Sina Bahram.

package Smuggle::Bind;
use strict; use warnings; use 5.008;

use CGI qw/:standard *table/;
use Data::Dumper;
use HTTP::Request::Common qw(POST);
use List::Util qw(min);
use LWP::UserAgent;

# Receives our POST requests.
my $URL = 'http://rna.williams.edu/cgi-bin/processNAR.cgi';

sub download {
    # Return the HTTP response from our POST request.
    my ($mrna, $tail) = @_;

    my $ua = LWP::UserAgent->new;
    my $req = POST $URL, [ whichDemo => '2',
                           string1 => $mrna,
                           string2 => $tail,
                         ];
    return $ua->request($req)->content;
}

sub scrape {
    # Takes the output from the POST request to BINDIGO and turns it
    # into a white-space delimited array of numbers. Raises a runtimer
    # error if unable to find the matrix.

    my ($content) = @_;
    $content =~ s/&nbsp;/ /g;

    # Match everything up to the `time` statistics.
    unless ($content =~ /(<br>\d.*)[ ]{9}.*real/s) {
        die('Smuggle::Bind: Unable to scrape matrix from BINDIGO');
    }

    $content = $1;
    $content =~ s/<br>/\n/gs;
    $content = substr($content, 1);
    return $content;
}

sub bind {
    # Returns an array of minimum energies along the diagonal of the
    # BINDIGO matrix. Raises a runtime error if mRNA or tail is
    # invalid or if the BINDIGO website is in an unexpected format.

    my ($tail, $mrna) = @_;

    # Some preconditions.
    if ($tail =~ /([^aucg])/) {
        die("Smuggle::Bind: Tail contains non-aucg letter $1");
    }

    if ($mrna =~ /([^aucg])/) {
        die("Smuggle:Bind: mRNA contains non-aucg letter $1");
    }

    if (length($mrna) < length($tail)) {
        die('Smuggle::Bind: mRNA is shorter than the tail');
    }

    # Here we go.
    my $n_tail = length($tail);
    my $len    = length($mrna) - $n_tail;

    return map {
        my $index    = $_;
        my $fragment = substr($mrna, $index, $n_tail);
        my $content  = download($fragment, $tail);
        my $matrix   = scrape($content);

        # A quick split and map to turn the list of numbers into a
        # matrix where a newline marks the start of each row.
        open(my $matrix_pipe, "+<", \$matrix)
          or die('Smuggle::Bind: Unable to open matrix pipe');
        my @array = map {[split /\s+/]} <$matrix_pipe>;
        close($matrix_pipe);

        # Sequences will conform to the minimal energy position, and
        # we only care about alignment of equal lengths. Remember that
        # $array[i][j] is the best way to align prefixes (1 .. i) with
        # (1 .. j). Therefore we only look at the minimum of the
        # matrix diagonal. Thanks, Vivek!
        my @diag = map {$array[$_][$_]} (1 .. $n_tail);
        my $min  = min @diag;

        # Works for Starmer, works for us.
        $min > 0 ? 0 : $min;
    } (0 .. $len);
}

sub bind_str {
    # Same as pull() except that it returns a MATLAB-friendly string
    # of floats instead of a Perl array.
    return join(' ', Smuggle::Bind::bind(@_));
}

1;