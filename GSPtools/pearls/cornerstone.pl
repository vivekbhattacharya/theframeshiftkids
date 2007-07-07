use strict; use warnings;
package Cornerstone;
use Smooth;

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

# The entire genome in a dainty little variable.
our $o;
sub extract {
    my ($start, $end) = @_;
    substr($o, $start, $end - $start + 1);
}

if (__FILE__ eq $0) {
    Smooth::helpcheck();
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
        print $handle Smooth::reverse_complement $_;
    } @sequences;
    
    print 'Done', $/;
}

1;
__END__


=head1 NAME

cornerstone.pl (web-enabled): Read a Genbank search file
and extract gene sequences from a genome file.

=head1 SYNOPSIS

=over

=item Script

=over 20

=item B<cornerstone.pl>

[B<--help>]

=item B<cornerstone.pl>

I<work folder> I<search results> I<genome file>

=back

=item Module

    use Cornerstone;
    Cornerstone::find_bounds $locations_file,
        sub { push @sequences, Cornerstone::extract(shift, shift) }
        sub { push @names, shift }

=back

=head1 DESCRIPTION

=over

=item $Cornerstone::o

This contains the entire genome as a string, all numbers
and cruft removed.

=item extract NUM, NUM

This returns the slice of $Cornerstone::o from the start
(first NUM) to the end (second NUM).

=item find_bounds STRING, SUB, SUB

This yields start and end boundaries of protein sequences
to the first SUB and the name of those proteins sequences to
the second SUB, both in the same order, after reading those
boundaries from the search results file stored in the path
STRING.

=end

=head1 OPTIONS AND ARGUMENTS

=over

=item I<work folder>

This folder becomes the base directory for specifying the
path to the genome file and the search results. In addition,
a new folder named C<proteins> will be created that stores
all the gene sequences.

=item I<search results>

=over

=item 1

Do a gene search for a string such as "NC_000913 heat shock."

=item 2

Show: 500. Send to: Text file.

=item 3

Save the text file as your search results. Cornerstone will
parse it. If you have more than 500 results, run cornerstone
on each one separately or concatenate the results (untested).

=back

=item I<genome file>

This comes from a Genbank Genome search. It should be a huge
list of nucleotides with or without numbers starting each line.

=back

=head1 DESCRIPTION

Cornerstone will look at all the proteins specified in the
search results and pull those sequences from the genome file.
Each gene sequence will be saved in a text file with filename
C<$protein.txt> in the C<proteins> folder mentioned above.

=head1 SEE ALSO

Genbank L<http://www.ncbi.nlm.nih.gov/>

=cut