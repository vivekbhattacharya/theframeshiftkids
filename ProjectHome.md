Click below for _good times_ by Hao Lian (shadytrees), Vivek Bhattacharya (nemosupport), and Daniel Vitek (drvitek), second place winners in the national Siemens Competition. You can download our presentation at the right. You can read the paper describing our model in lustrous detail with full references to biological literature as well as a cursory synthesis of physics, biology, mathematics, and angst.

[![](http://lh5.google.com/shadytrees/Rrdz5hPfOTI/AAAAAAAAANY/RNtTnz1TFOs/thrush.jpg)](http://theframeshiftkids.googlecode.com/files/GSPtools-4.669.zip)

[GSPtools is documented within the work-in-progress manual.](Manual.md)

## News ##
See [Changelog](Changelog.md) for details.
  * April 10, 2008: GSPtools 4.669 is out, and "fast" can barely describe it.
  * You can now view the mathematical documentation behind the model in the "Siemens Paper.pdf" file to your right.
  * November 4: GSPtools 4.66 is out. You will need to download the new Travel2.mat if you want to use the new values developed by our genetic algorithm, affectionately named God 2.0. If not, rename Travel.mat to Travel2.mat. By the way, our project won regionals at the Siemens Competition.
  * August 23: GSPtools 4.6 is out. Changelog updates are forthcoming. The largest change is that we've reverted to 3.141's non-sensitive model, effectively reversing [the previous sensitivity changes](BecomingHarsher.md). However, we are keeping with the 4.x versioning due to the major structural changes engendered in 4.
  * July 15: GSPtools 4 released with major structural and algorithmic changes. Consult the Changelog. This should be perceived as an _entirely new model_ that is not compatible with the old model. For example, prfB has a yield of 10% with this model right now. (We are investigating why.)
  * July 1: GSPtools 3.141 released with huge performance increases. Remember to replace `Codons.mat` and `TAV.mat` with `Travel.mat`. And always delete your old GSPtools folder first.

## Summary ##
GSPtools (genetic signal processing tools) is a toolbox of Matlab and Perl code originally developed by [Dr. Lalit Ponnala](http://sites.google.com/site/lalitp/) and then rewritten by us. It takes mRNA nucleotide sequences (FASTA or not) and produces a displacement plot and a polar plot of the phase of the free energy signal. That's not all: Also included is a genetic algorithm that calibrates tRNA availability values, a greedy algorithm to compute an optimal synonymous mRNA sequence, sensitivity plots to test model parameters, a Perl program that pulls Ecogene genes, and another Perl program that pulls genes from a a genome file. And we haven't even gotten to the completely rewritten, modular free energy Perl module Kidnap based on [Joshua Starmer's free2bind project](http://sourceforge.net/projects/free2bind) or the metrics we've developed to calculate protein yield from an mRNA sequence.

New free energy calculations come from the [BINDIGO](http://rna.williams.edu/tools/BINDIGO.html) algorithm written by Nathaniel Hodas and [Professor Daniel Aalberts](http://panic.williams.edu/).  As per his request, code to run the BINDIGO algorithm is not available from this webpage.  Please contact Professor Aalberts to request permission to obtain the code related to BINDIGO as part of our program.

## Quick Start ##

  * Download and install [Strawberry Perl](http://strawberryperl.com/).
  * Choose a directory (we'll call it `C:\Work` here.)
  * Download and unzip the latest version of our code (see `Featured Downloads`) to `C:\Work\GSPtools`.
  * Copy the below prfB sequence to a file (we'll call it `prfb.fasta` here) in `C:\Work\`.
  * Download `Travel.mat` to `C:\Work\`.
  * Open Matlab. Go to `File -> Add Paths -> Add Subdirectory -> Add C:\Work.`

### Having fun ###
  * For one iteration, type `unity('prfb.txt', 25)` into Matlab.
    * 25 because there's an intended frameshift at codon 25.
  * For a continuous loop to measure average yield (Law of Large Numbers, wink wink), type `megaunity('prfb.txt', 25)`.
  * For Internet fun, type `unity('http://shadytrees.pastebin.ca/raw/538494', 25)`.
    * Or `megaunity('http://shadytrees.pastebin.ca/raw/538494', 25)`

Of course, the model works for other E. coli sequences too. We have tested the entire genome, ribosomal proteins in _E. coli_, Weiss genes, bovine growth hormone, and our artificial frameshifter.

```
> prfB from E. ravioli
agaaaucagaccauguuugaaauuaauccgguaaauaaucgcauucaggaccucacggaa
cgcuccgacguucuuaggggguaucuuugacuacgacgccaagaaagagcgucuggaaga
aguaaacgccgagcuggaacagccggaugucuggaacgaacccgaacgcgcacaggcgcu
ggguaaagagcguuccucccucgaagccguugucgacacccucgaccaaaugaaacaggg
gcuggaagauguuucuggucugcuggaacuggcuguagaagcugacgacgaagaaaccuu
uaacgaagccguugcugaacucgacgcccuggaagaaaaacuggcgcagcuugaguuccg
ccguauguucucuggcgaauaugacagcgccgacugcuaccucgauauucaggcgggguc
uggcgguacggaagcacaggacugggcgagcaugcuugagcguauguaucugcgcugggc
agaaucgcgugguuucaaaacugaaaucaucgaagagucggaaggugaaguggcggguau
uaaauccgugacgaucaaaaucuccggcgauuacgcuuacggcuggcugcguacagaaac
cggcguucaccgccuggugcguaaaagcccguuugacuccggcggucgucgccacacguc
guucagcuccgcguuuguuuauccggaaguugaugaugauauugauaucgaaaucaaccc
ggcggaucugcgcauugacguuuaucgcacguccggcgcgggcggucagcacguuaaccg
uaccgaaucugcggugcguauuacccacaucccgaccgggaucgugacccagugccagaa
cgaccguucccagcacaagaacaaagaucaggccaugaagcagaugaaagcgaagcuuua
ugaacuggagaugcagaagaaaaaugccgagaaacaggcgauggaagauaacaaauccga
caucggcuggggcagccagauucguucuuauguccuugaugacucccgcauuaaagaucu
gcgcaccgggguagaaacccgcaacacgcaggccgugcuggacggcagccuggaucaauu
uaucgaagcaaguuugaaagcaggguuauga
```