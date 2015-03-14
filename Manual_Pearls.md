Besides free energies with `scan_brightly.pl`, you can do a smorgasbord of other bioinformatics-y things with the pearls including, but not limited to,

  * Have a `scan_brightly.pl` sleepover (not documented here)
  * Most importantly, download all the genes from another genome
  * Create a `Travel.mat` for _any_ organism
  * Convert a string of amino acids, written in the one-letter abbreviation form to an mRNA sequence
  * Test the pearls

Besides the in-depth descriptions below, all Perl programs in BWFtools should have built-in POD documentation. That is, if you type the program at the command-line without any arguments (e.g. just `dysentery.pl`), you should get a steamy short story about an illicit relationship between panthers and mountain lions and/or helpful documentation.

## Downloading genes ##
```
> cd "/path/to/BWFtools/pearls"
> echo NC_000913| perl pull_genome.pl my_genomes
> type list_of_genes.txt | perl pull_genome.pl my_genomes
> genbankest.pl my_genomes
```

Ingredients required: You'll need Internet access, specifically access to [Genbank's website.](http://www.ncbi.nlm.nih.gov/) You'll also need the accession number for your organism. For example, _E. coli_'s K12 strain corresponds to NC\_000913. You can do a quick search of Genbank's website and hope that the weird algorithm they use proffers you right one. An alternative is to search "Genbank _organism name_" on Google, which yields better results. Let's assume your accession number is A1821.

To pull the genome, open a command prompt to the `pearls` folder and type `echo A1821|perl pull_genome.pl my_genomes`. This will save a `A1821.txt` file to the `my_genomes` folder. Now type `genbankest.pl my_genomes`. Genbankest downloads the [Nucleotide file](http://www.ncbi.nlm.nih.gov/sites/entrez?db=nucleotide) for `A1821.txt` and saves it to `A1821-genbankest.txt` for internal use by `genbanker.pl`. Afterward, you should see all the available genes in the `my_genomes/A1821` directory.

To pull multiple genomes, open a text file and type each one on a new line. Save it to the `pearls` folder. Let's assume you saved it as `my_genomes.txt`. At the command line, type `type my_genomes.txt | perl pull_genome.pl my_genomes`. `genbankest.pl my_genomes` works the same. The genes for each genome will be saved to a different director in `my_genomes`.

Genbank's Nucleotide files are, I surmise, designed by a somewhat-intelligent marshmallow; they are notoriously difficult to parse. `genbanker.pl`, which is the subsystem for `genbankest.pl`, is friendly and verbose in reporting when it cannot pull a certain gene at which point you can send us (and by us, we mean shadytrees) an email with the accession number and gene. (The flip side is that `genbanker.pl` might not catch the error and _everything_ will go haywire.)

Nucleotide files also sometimes give loci the same gene name. For example, there might be two loci (say, 100..200 and 200..300) with the gene name "H" because they both produce "H-protein," and somebody thought it'd be super-clever to have the end-user work through this edge case. We deal with it by saving the first one to `H.txt` and the rest to `H-again-x.txt` where "x" is some random number chosen to prevent conflicts.

## Creating wait times ##
```
> codonusagetable2freqs.pl http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=602 > /path/typhimurium-freqs.txt
M> Travel = freqs2travel('/path/typhimurium-freqs.txt');
M> save('/path/Travel.mat', 'Travel');
```

Beached Whale Frameshift tools are built to run _E. coli_ mRNA sequences, whose codons have wait times specific to _E. coli_. In addition to changing the species angle and the C1 brownie constant in `displacement.m`, you must also change the wait times (living within `Travel2.mat`) once you work on another organism. For example, _Salmonella typhimurium_. Although [codon usage information](http://www.kazusa.or.jp/codon/) is freely available online, wait times are not. However, `freqs2travel.m` converts codon usage frequencies to a usable `Travel` structure.

Ingredients: URL to a codon usage table. For example, [the URL for S. typhimurium](http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=602). HTML is not a suitable format for `freqs2travel.m`, but you can resolve this mini-crisis with

```
codonusagetable2freqs.pl http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=602 > /path/typhimurium-freqs.txt
```

Now open Matlab and type `Travel = freqs2travel('/path/typhimurium-freqs.txt')`. You now have the structure you need stored in a Matlab variable. To save that variable as a file, type `save('/path/Travel.mat', 'Travel')`. Note that you type Travel as a string, not as a variable identifier.

Now, if you perform this with _E. coli_ K12, the results will not match with `Travel2.mat`. This is because we optimized our wait times for prfB using the God 2.0 algorithm on the shaky basis that evolution would produce a ratio of tRNA cells in the cell optimal for translating the most often used proteins; for lack of time to find and run the most used proteins (and a lack of research in this area), we used prfB as the gold standard.

## Going from amino acids to mRNA and back ##

Sometimes, for various socioeconomicpolitic reasons, you aren't given a mRNA sequence. Before you start the revolution, realize that BWFtools has stepped in to help your thorny problems of nitrogen, phosphorous, and _el amor_. The program in question is `dysentery.pl`. Because it's a straightforward program and its documentation is ship-shape, this is a walk-through of one of its **three** possible functions.

Ingredients: a plain-text file of a string of amino acids. You should use the standard one-letter abbreviations them, which is available not only on [Wikipedia](http://en.wikipedia.org/w/index.php?title=Proteinogenic_amino_acid&oldid=224377445#Chemical_properties) but also in `Smooth.pm`. For example, say you just hand-copied the amino acid sequence for human insulin from a high-school biology textbook into `insulin.txt`. To convert, you'd run `dysentery.pl --run insulin.txt auauauauauau` where `auauauauauau` is an imaginary 12-nucleotide leader sequence you'd need to run it in _E. coli_.

Caveat: If you'll hark back to said biology course for a moment, amino acids have more than one codon correspondence. `dysentery.pl` will simply randomly choose one. Therefore, you will receive different mRNA sequences for each run. Feature to make up for the caveat: Dysentery is intelligent about handling URLs like files, meaning instead of typing "insulin.txt" above you could have typed "http://example.org/insulin.txt" and it still would've worked. For more information about the other two modes and working example URLs, consult the Dysentery documentation.