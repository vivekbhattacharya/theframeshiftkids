## Defines some substr calls used in force_bind. See
## force_bind() for more details.
##
## Note: fivethree is threefive for a y doublet and
## similarly ->{context} is switched.
use strict; use warnings;
package Strand;
# Asterisk needed to workaround substr's handling
# of negative numbers. That is, we want $i-1 to
# still work even when it doesn't produce s55.
sub new {
	my ($class, $seq, $length) = @_;
	my $self = {seq => '*'.$seq, length => $length};
	
    #helper($self, $i);
	bless $self, $class;
}

sub update {
    my ($self, $i) = @_;
    my ($s55, $s5, $s3, $s33) = split '', substr($self->{seq}, $i, 3);
	$self->{fivethree} = [$s5, $s3];
	$self->{context} = ($i > 0 && $i < $self->{length}-2) ? 
		[$s55, $s33] : ['', ''];
}

sub all { @{shift->{fivethree}} }
sub context { @{shift->{context}} }

1;