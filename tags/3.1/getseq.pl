# Extracts ONLY the character sequence from manually-created file
use Smooth qw(getseq);
print Smooth::getseq(shift @ARGV);