# Extracts ONLY the character sequence from manually-created file
use Smooth;
my $seq = Smooth::getseq($ARGV[0]);
print $seq;

open(my $handle, qq/>$ARGV[1]/);
print $handle '>sequence', $/, $seq;