#!/usr/bin/perl

##========================================================================
##
## Author: Joshua Starmer <jdstarme@unity.ncsu.edu>, 2004
## 
## Copyright (C) 2004, Joshua Starmer
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
## 
##========================================================================

use warnings;
use strict;

use POSIX qw(tmpnam); # this generates filenames for temporary files
use Fcntl; # this defines the constants used in sysopen
use Getopt::Std; 
use FreeAlign;


##========================================================================
##
## DEFINE GLOBAL CONSTANTS
##
##========================================================================

use constant TRUE => 1;
use constant FALSE => 0;

use constant SINGLE_COLUMN => 1;
use constant MULTI_COLUMN => 2;
use constant MULTI_ROW => 3;
use constant GRAPH_MODE => 4;


##========================================================================
##
## DEFINE GLOBAL VARIABLES
##
##========================================================================

my $gDebugMode = FALSE;

my $gQuietMode = FALSE;

my $gEnergyOnly = FALSE;

my $gOutputMode = SINGLE_COLUMN;

my $gLeaderLength = 0;
my $gTrailerLength = 0;
my $gNumBases;

my $gDoubletID = "XIA_MATHEWS";
my $gLoopID = "JAEGER";

my $gTerminalBulgeRule = "ALLOWED";
my $gInternalBulgeRule = "ALLOWED";
my $gLoopRule = "ALLOWED";

my $gTemperature = 37;

my $gRnaSeq;
my $gSeqFile;
my $gIndexData;



##========================================================================
##
## DEFAULT HELP MESSAGE
##
##========================================================================
my $sHelp = <<END_OF_HELP;

Usage: free_scan.pl [options] rna_seq sequence.fasta
    
    free_scan takes an RNA binding sequence (like the 3 prime tail on the 16s 
    rRNA) and sees how well it will bind to the sequence (DNA or RNA, if DNA
    it is first converted to RNA) at each point in a range specified in 
    index_data.
    
    All of the sequence contained in sequence.fasta will be examined.

    NOTES:
    The sequence data contained in "sequence.fasta" is assumed to be
    in 5'->3' orientation.

    "rna_seq" is assumed to be written in 3'->5' orientation.

    example:
    shell> ./free_scan.pl -p XIA_MATHEWS auuacuagg human_dna.fasta human_gene_indices > delta_G_values.txt

    Commonly used strings for "rna_seq" are 3 prime 16/18 S tails.
    Here are several listed as reference, written 3prime->5prime:

    E. coli: auuccuccacuag (13 bases)
    Z. mays: guuacuag (8 bases)
    Arabidopsis thaliana: guuacuag (8 bases)
    Xenopus laevis: auuacuag (8 bases)
    Mus musculus: auuacuag (8 bases)
    D. melanogaster: auuacuag (8 bases)
    H. sapiens: auuacuag (8 bases)

Options:

-h
   print this help message.

-d
   Debug Mode.  Prints out all sorts of debugging information...

-c
   Print the output in multiple columns separated by tab characters.  The
   default is to print the output as a single column

-r
   Print the output in multiple rows.

-e
   Just print energies.  Refrain from printing comments in the output.

-g
   Print the results in two columns, with relative spacing (relative to the
   0 position) in the first column and Delta G values in the second column.
   This option assumes that the rna_seq is the 3\' 16S rRNA tail.

-p FREIER|SANTALUCIA|XIA_MATHEWS
   Allows you to determine which set of parameters are used to simulate binding
   between the to strands of nucleotides.  The default value is XIA_MATEHWS.
   FREIER - RNA binding parameters from 1986
   SANTALUCIA - DNA binding parameters from 1998
   XIA_MATHEWS - RNA binding parameters from 1998

-t temperature
   The temperature, in Celsius, at which the binding is supposed to take place.
   The default value is 37.

-q
   Quiet Mode.  Refrain from printing out progress messages.  Overrides Debug
   Mode if it is also selected.
END_OF_HELP


##========================================================================
##
## SUBROUTINE print_headers(): print output headers
##
##========================================================================
sub print_headers {
    my ($delimiter, $f_handle) = @_;

    if (!defined($delimiter)) {
	$delimiter = "\n";
    }

    if (!defined($f_handle)) {
	$f_handle = *STDOUT;
    }

    print $f_handle "\# Leader Length: ".$gLeaderLength.$delimiter;
    if ($gNumBases) {
	print $f_handle "\# Num Bases: ".$gNumBases.$delimiter;
    }
	print $f_handle "\# Force Bind: ON".$delimiter;
    print $f_handle "\# Doublet Binding Parameters: ".$gDoubletID.$delimiter;
    print $f_handle "\# Internal Loop Parameters: ".$gLoopID.$delimiter;
    print $f_handle "\# Terminal Bulges: ".$gTerminalBulgeRule.$delimiter;
    print $f_handle "\# Internal Bulges: ".$gInternalBulgeRule.$delimiter;
    print $f_handle "\# Internal Loops: ".$gLoopRule.$delimiter;
    print $f_handle "\# Binding Temperature: ".$gTemperature." C".$delimiter;
    print $f_handle "\# Binding Sequence: 3'-".$gRnaSeq."-5'".$delimiter;
    print $f_handle "\# mRNA Sequence File: ".$gSeqFile.$delimiter;	
}


## "main"
{
    
    # check to see if anything is coming in to the program (via a pipe or
    # on the command line).  If there is no input at all, print out the help
    if (- STDIN && !@ARGV) {
        print STDERR $sHelp;
        exit 1;
    }
    
    # process the command line arguments
    my %opts; # a hash table to store file names passed in as agruments
    getopts('ehdqcL:T:n:fp:t:Ors:g', \%opts);
     
    if ($opts{'h'}) { # print help and exit
        print STDERR $sHelp;
        exit 1;
    }

    if ($opts{'d'}) { # turn on debugging mode!
	$gDebugMode = TRUE;
    }

    if ($opts{'q'}) { # turn on quiet mode!
	$gQuietMode = TRUE;
	$gDebugMode = FALSE;
    }

    if ($opts{'e'}) { # don't print comments in output file
	$gEnergyOnly = TRUE;
    }

    if ($opts{'c'}) { # print the output in multiple columns instead of just 1
	$gOutputMode = MULTI_COLUMN;
    }
    
    if ($opts{'r'}) { # print the output one row per row in index_data
	$gOutputMode = MULTI_ROW;
    }

    if ($opts{'g'}) { # print the output in two columns (RS and Delta G)
	$gOutputMode = GRAPH_MODE;
    }

    if ($opts{'n'}) {
	$gNumBases = $opts{'n'};
    }

    # now that we know we are not just going to exit the program, create
    # an FreeAlign object...
    my $align_obj = new FreeAlign();

    if (defined($opts{'p'})) {
	$gDoubletID = $opts{'p'};
    }

    if ($gDoubletID eq "FREIER") {
	$align_obj->select_parameters($align_obj->FREIER);
    } elsif ($gDoubletID eq "SANTALUCIA") {
	$align_obj->select_parameters($align_obj->SANTALUCIA);
    } elsif ($gDoubletID eq "XIA_MATHEWS") {
	$align_obj->select_parameters($align_obj->XIA_MATHEWS);
    } else {
	print STDERR "ERROR: gDoubletID : ".$gDoubletID." undefined\n";
	exit;
    }    

    if (defined($opts{'t'})) {
        $gTemperature = $opts{'t'};
        $align_obj->{Temp} = $gTemperature;
    }

    # turn off loops and bulges...
    $align_obj->{InternalBulge} = FALSE;
    $align_obj->{Loop} = FALSE;

    # get the rna binding sequence
    $gRnaSeq = shift(@ARGV);
    if (!defined($gRnaSeq)) {
	print STDERR "\nERROR: you must specify an rna binding sequence\n\n";
        print STDERR $sHelp;
	exit;
    }
    $gRnaSeq = $align_obj->dna2rna($gRnaSeq);
    $gRnaSeq = uc($gRnaSeq);
    

    # get the sequence file name
    $gSeqFile = shift(@ARGV);
    if (!defined($gSeqFile)) {
	print STDERR "\nERROR: you must specify a sequence file\n\n";
        print STDERR $sHelp;
	exit;
    }

    my $path = shift(@ARGV);

    if (!$gQuietMode) {
	print STDERR "Extracting sequence data from $gSeqFile.\n";
    }
    
    # print out header information for the output...
    if (!$gEnergyOnly  && !(defined($path))) {
	if ($gOutputMode != MULTI_ROW) {
	    print_headers();
	    if ($gOutputMode == GRAPH_MODE) {
		print "# RS\tDelta G\t#Absolute Pos.\tBase\n";
	    }
	} else {
	    print_headers("\t");
	}
    }

    my ($full_seq) =  $align_obj->load_fasta($gSeqFile);

    my $first_pass = TRUE;
    my @output_matrix;
    my $output_cols=0;
    my $counter = 1;
    my @index_data = ("ALL");
    foreach my $starts_and_stops (@index_data) {
	my $seq_data = $full_seq;
	chomp($seq_data);
	
	$seq_data = $align_obj->dna2rna($seq_data);
	$seq_data = uc($seq_data);
	
	if (!$gQuietMode) {
	    print STDERR "\n$counter: Starting to analyze sequence data: $starts_and_stops\n";
	    $counter++;
	}
	
	my $rna_binding_length = length($gRnaSeq);
	my $rs_offset = $rna_binding_length - 5;
	my $seq_length = length($seq_data);


	if ($first_pass && ($gOutputMode == MULTI_ROW)) {
	    $first_pass = FALSE;

	    for (my $i=0; $i<=($seq_length - $rna_binding_length); $i++) {
    
            print ($i-$gLeaderLength);
    
            if ($i != ($seq_length - $rna_binding_length)) {
                print "\t";
            } else { print $/ }
	    }
	}
	
	my $num_loops = $seq_length - $rna_binding_length;
	my $time_check = $num_loops / 10;
	my $timer = 0;
	for (my $i=0; $i<=($seq_length - $rna_binding_length); $i++) {
	    
	    if (!$gQuietMode) {
		$timer++;
		
		if ($timer > $time_check) {
		    print STDERR ".";
		    $timer = 0;
		}
	    }
	    
	    my $seq = substr($seq_data, $i, $rna_binding_length);
	    
	    my $first_char = substr($seq, 0, 1);

		$align_obj->force_bind($seq, $gRnaSeq);
	    # get the score for the best binding...
	    my $free_energy = $align_obj->total_free_energy();
	    
	    # free energy values greater than zero represent binding that would
	    # not take place (without added energy) and thus, should be set to
	    # zero (the same as if no binding had happened at all).
		$free_energy = 0 if $free_energy > 0;
	    
	    my $out_string = $free_energy;
	    
	    if (!$gEnergyOnly && !($gOutputMode == MULTI_ROW)) {
		$out_string 
		    = $out_string."\t# ".$first_char." ".($i-$gLeaderLength);
	    }

	    if ($gOutputMode == SINGLE_COLUMN) {
		if (defined($path)) {
		    print OUT_FILE $out_string."\n";
		} else {
		    print $out_string."\n";
		}
	    } elsif ($gOutputMode == MULTI_COLUMN) {
		$output_matrix[$output_cols][($i+1)] = $out_string;
	    } elsif ($gOutputMode == MULTI_ROW) {
		if (defined($path)) {
		    print OUT_FILE "\t".$out_string;
		} else {
		    print "\t".$out_string;
		}
	    } elsif ($gOutputMode == GRAPH_MODE) {
		if (defined($path)) {
		    print OUT_FILE "".($i-$gLeaderLength+$rs_offset)."\t".$free_energy."\t#".($i-$gLeaderLength)."\t".$first_char."\n";
		} else {
		    print "".($i-$gLeaderLength+$rs_offset)."\t".$free_energy."\t#".($i-$gLeaderLength)."\t".$first_char."\n";
		}
	    }
	    

	}

	if ($gOutputMode == SINGLE_COLUMN) {
	    if (defined($path)) {
		print OUT_FILE "\n";
	    } else {
		print "\n";
	    }
	} elsif ($gOutputMode == MULTI_COLUMN) {
	    $output_cols++;
	} elsif ($gOutputMode == MULTI_ROW) {
	    print "\n";
	}
    }

    if (($gOutputMode == SINGLE_COLUMN) || ($gOutputMode == MULTI_ROW)
	|| ($gOutputMode == GRAPH_MODE)) {
	if (!$gQuietMode) {
	    print STDERR "\nDone!\n";
	}
	exit 1; # we've already printed things out so we are done!
    }


    if (!$gQuietMode) {
	print STDERR "\nPrinting Output!\n";
    }

    # print everything out...
    # The way the array is indexed is a little odd here, you can get
    # details in "man perldsc" (man page for perl data structures)
    my $max_value = 0;
    my $max_index = 0;
    for my $h (0 .. $output_cols) {
	if ($#{$output_matrix[$h]} > $max_value) {
	    $max_value = $#{$output_matrix[$h]};
	    $max_index = $h;
	}
    }

    for my $i (0 .. $#{$output_matrix[$max_index]}) {
	for my $h (0 .. $output_cols) {
	    if (defined($output_matrix[$h][$i])) {
		print $output_matrix[$h][$i]."\t";
	    } else {
		if (!$gEnergyOnly) {
		    print "\t\t";
		} else {
		    print "\t";
		}
	    }
	}
	print "\n";
    }

    if (!$gQuietMode) {
	print STDERR "\nDone!\n";
    }

    exit 1;
}
