# Extracts ONLY the character sequence from manually-created file
use strict;

sub getseq {
   my $file = shift; {
      open(my $handle, $file) or die "getseq cannot open \"$file\" for reading\n\n";
      for(join '', <$handle>) { s/[\s0-9]//g and return $_; }
   }
}

print getseq(shift @ARGV) unless caller;