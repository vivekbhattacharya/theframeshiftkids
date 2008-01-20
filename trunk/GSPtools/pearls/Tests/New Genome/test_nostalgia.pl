use strict;
use warnings;
use lib '../..';
use Smooth qw(getseq);

my ($old_dir, $new_dir) = @ARGV;

chdir($old_dir);
my @old = glob("*.txt");
chdir('..');
while (1) {
    my $file_old = $old[rand @old];
    my $file_new = "$new_dir/$file_old";
    $file_old = "$old_dir/$file_old";

    next if (!-f $file_new);

    my $seq_old = Smooth::getseq $file_old;
    my $seq_new = Smooth::getseq $file_new;

    if ($seq_old eq $seq_new) {
        # print "$file_old: ok\n";
    }
    else {
        print "$file_old: FAIL\n";
    }
}
