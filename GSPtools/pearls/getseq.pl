$help = <<END;
NAME
    getseq.pl
    (web-enabled)

USAGE
    getseq.pl [path to sequence] [path to fasta file]
    
    getseq.pl retrieves the sequence and prints it.
    In addition, it saves the sequence in a FASTA
    format suitable for free2bind.
END
if (!@ARGV or $ARGV[0] eq '--help') { print $help; exit }
   
# Extracts ONLY the character sequence from manually-created file
use Smooth;
my $seq = Smooth::getseq($ARGV[0]);
print $seq;

# Write sequence to a fasta format.
open(my $handle, qq/>$ARGV[1]/);
print $handle '>sequence', $/, $seq;