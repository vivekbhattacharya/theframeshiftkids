# Usage: `cornerstone.pl folder locations genome`
# Output: a 'proteins' folder of files with protein sequences
###  folder: work folder and possibly the location of
###      `locations` and genome; it is the
###      place where cornerstone will place the "proteins" folder
### locations: file containing the list of sequence locations
###      generated by Genbank's search "To Text" function;
###      by default, cornerstone will search in the folder you specified
###      for `folder`, which you can override with a full path
### genome: haystack to search (e.g. the _E. coli_ K12 genome);
###      like `locations`, cornerstone will first search in
###      the work folder you specified

use strict;
use warnings;
use Smooth;
package Cornerstone;

# This is a generator. It yields beginning and ending
# indices first, then the protein names.
sub find_bounds {
    my ($file, $push_bounds, $push_name) = @_;
    Smooth::webopen $file, sub {
        if ($_ =~ /(\d*)\.\.(\d*)/) {
            my ($start, $end) = ($1, $2);            
            # Quirk: Genbank sometimes reverses the indices
            ($start, $end) = ($end, $start) if ($start > $end);
            # Leader sequence of 12 before `tac/uag` starts
            &$push_bounds($start, $end + 12);
        } elsif ($_ =~ /^\d*?: (.*)$/) {
            my $hair = $1; # Vivek's hairy
            $hair =~ s/\s//g;
            &$push_name($hair);
        }
    }
}

# Cf. page 43-6 of _Programming Perl for Bioinformatics_
# by James D. Tisdall
sub reverse_complement {
    for(shift) { tr/atcg/uagc/; return scalar reverse; }
}

# The entire genome in a dainty little variable.
our $o;
sub extract {
    my ($start, $end) = @_;
    substr($o, $start, $end - $start + 1);
}

our $help = <<END;
NAME
    cornerstone.pl
    (web enabled)

USAGE
    perl cornerstone.pl [work directory]
        [list of locations (file name)]
        [genome (file name)]
    
    Retrieves the gene sequences specified in a Genbank
    search results (text file format) from the genome
    given. Stores each gene sequence in a text file under
    a 'proteins' directory in the working directory. Each
    text file is named after the gene sequence name e.g.
    prfB.txt.
END
if (__FILE__ eq $0) {
    if (!@ARGV or $ARGV[0] eq '--help') { print $help; return }
    my ($dir, $locations, $genome) = @ARGV;
    chdir($dir);
    
    # Slurp the E. Coli genome as a huge text file
    # for quick access.
    $o = Smooth::getseq($genome);
    
    my (@sequences, @names);
    # The first anonymous function pushes yielded bounds
    # into the @sequences array after extracting the
    # sequences with `extract`.
    # The second anonymous function pushes it
    # without any processing.
    find_bounds $locations,
        sub { push @sequences, &extract; },
        sub { push @names, shift; };
    
    # Output: store to file with special filenames
    mkdir('proteins'); chdir('proteins');
    map {
        my $name = shift @names;
        open(my $handle, ">$name.txt") or
            die("Cornerstone cannot open $name for writing");
        print $handle reverse_complement $_;
    } @sequences;
    
    print 'Done', $/;
}

1;