## Note: fivethree is threefive for a y doublet and
## similarly ->{context} is switched.
use strict; use warnings;
package Strand;

# Asterisk needed to workaround substr's handling
# of negative numbers. That is, we want $i-1 to
# still work even when it doesn't produce s55 and we
# want $i+3 to yield a non-null value for $s33 even
# when none exists.
sub new {
	my ($class, $seq, $length) = @_;
	my $self = {};
    $self->{seq} = [split q//, qq/*$seq*/];
    $self->{length} = $length;
	
	bless $self, $class;
}

sub update {
    my ($self, $i) = @_;
    my ($s55, $s5, $s3, $s33) = @{$self->{seq}}[$i .. $i+3];
	
    $self->{all} = [$s5, $s3];
    $self->{context} = $i > 0 && $i < $self->{length}-2 ? 
		[$s55, $s33] : ['', ''];
}

sub all { @{shift->{all}} }
sub context { @{shift->{context}} }

1;