
#!/usr/bin/perl
# Extracts ONLY the character sequence from manually-created file

use strict;
my $file = shift @ARGV; {
        open(my $handle, $file) or die "Cannot open file \"$file\" for reading!\n\n";
        $_ = join '', <$handle> and s/[\s0-9]//g and print;
}