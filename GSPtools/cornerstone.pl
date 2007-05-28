use strict;
package Cornerstone;

sub find_bounds {
	my ($file, $push_bounds, $push_name) = @_;
	open(my $handle, $file) or die "Cannot open the file \"$file\" for reading!\n";
	while(<$handle>) {
		if ($_ =~ /(\d*)\.\.(\d*)/) {
			my ($start, $end) = ($1, $2);
			($start, $end) = ($end, $start) if ($start > $end);
			&$push_bounds($start, $end + 12);
		} elsif ($_ =~ /^\d*?: (.*)$/) {
			my $hair = $1; # Vivek's hairy
			$hair =~ s/\s//g;
			&$push_name($hair);
		}
	}
}

our $o; prep();
sub prep {
	open(my $handle, 'ecoli-k12.txt') or die "Cannot open Ecoli genome\n";
	while(<$handle>) {
		s/\s//g and /^(\d*)(.*)$/;
		$o .= $2;
	}
}

sub extract {
	my ($start, $end) = @_;
	return substr($o, $start, $end - $start + 1);
}

# Page 43-6 of _Programming Perl for Bioinformatics_ by James D. Tisdall
sub reverse_complement {
	$_ = shift and tr/atcg/uagc/ and return scalar reverse $_;
}

########## Main Logic ##########
our(@sequences, @names);
find_bounds(@ARGV, sub { push @sequences, &extract; },
	sub { push @names, shift; });

mkdir('proteens');
chdir('proteens');
map {
	my $name = shift @names; {
		open(my $kampf, ">$name.txt") or die("Cannot open $name for writing");
		print $kampf reverse_complement $_;
	}
} @sequences;