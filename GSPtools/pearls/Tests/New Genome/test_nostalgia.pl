use strict; use warnings;

use File::Basename;
use lib dirname(__FILE__) . '/../..';
use Smooth qw(getseq);

my ($old_dir, $new_dir) = @ARGV;

my @old = glob($old_dir . "/*.txt");
while (1) {
    my $file_old = basename $old[rand @old];
    my $file_new = "$new_dir/$file_old";
    $file_old = "$old_dir/$file_old";

    if (!-f $file_new) {
        print "I can't find $file_old in $new_dir\n";
        next;
    }

    my $seq_old = Smooth::getseq $file_old;
    my $seq_new = Smooth::getseq $file_new;

    if ($seq_old eq $seq_new) {
        # print "$file_old: ok\n";
    }
    else {
        print "$file_old: FAIL\n";
    }
}
