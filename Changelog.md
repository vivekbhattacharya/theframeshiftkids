## What's New in Version Four ##

### 4.6692 ###
  * It can now single-handedly [extract all genes from any genome.](Manual_Pearls.md)
  * More vectorization, especially `inst_energy.m` (`diff_vectors`), means more performance
  * Can convert [Codon Usage Database](http://www.kazusa.or.jp/codon/) entries into `Travel`s (loop times)

### 4.669 ([r643](https://code.google.com/p/theframeshiftkids/source/detail?r=643)) ###
  * We now compute `inst_energy` with actual vector subtraction rather than converting between phase-magnitude and complex numbers.
  * It's faster than ever thanks because Kidnap now caches free energy results.
  * We renamed `cumm_mag_phase` and `diff_vectors` to `cumm_energy` and `inst_energy` to reflect their actual use.
  * We removed stop detection due to its being obsolete.
  * Performance increases as we continue vectorization and after buffering `disp_shifts.m` and implementing quick-fail `Config.dire`.
  * `quasiforce` plots the wait times of a sequence. Don't ask me why.

### 4.66 ([r524](https://code.google.com/p/theframeshiftkids/source/detail?r=524)) ###
  * God 2.0: It's our new genetic algorithm that attempts to find the best TAV matrix given a yield criteria. In our case, we used it to find the optimal TAV values in an attempt to calibrate them according to bovine growth hormone biological data.
  * Performance increases across the board as we vectorize GSPtools.
  * Config.detect\_stops and Config.detect\_pauses controls the detection of stop codons and long wait times.
  * sensitivity.m creates sensitivity plots, mainly used for Weiss genes.
  * Travel2.mat replaces the old Travel.mat if you want the new values outputted from God 2.0.

### 4.6 ([r325](https://code.google.com/p/theframeshiftkids/source/detail?r=325)) ###
  * [The harsher yield model](BecomingHarsher.md) is now gone, but in its place comes [Deviation](Deviation.md). You can switch back to probability via `config.m`.
  * Extensive documentation for the `Kidnap` module ([Issue 3](https://code.google.com/p/theframeshiftkids/issues/detail?id=3))
  * Optional arguments for `megaunity` and `unity` ([Issue 1](https://code.google.com/p/theframeshiftkids/issues/detail?id=1)); for usage, see in-file documentation.
  * Model now detects stop codons and unusually long wait times (`displacement.m`) that we call ribosomal hyperpauses for short
  * `alopeciaunity.m` saves polar plots in batches, similar to what `opportunity.m` and `jejunity.m` do for overlay and error bar plots
    * `impunity.m` handles single files in addition to folders now. Try it out with `impunity('prfb.txt', {'uga,25'}, {}, 100)`.

### 4 ([r212](https://code.google.com/p/theframeshiftkids/source/detail?r=212)) ###
  * reloop and sofar gone, q.v. BecomingHarsher
  * The ThrushBaby algorithm

  * free2bind, refactored to Kidnap, heavily edited (from 9000 lines to 400)
    * Now extensible with other parameter modules and readable
    * SantaLucia module depreciated due to lack of key subroutines
    * free2bind.zip no longer a requirement for installation
  * Performance increase in `diff_vectors.m`

## What's New in Version Three ##

### 3.141 ([r156](https://code.google.com/p/theframeshiftkids/source/detail?r=156)) ###
  * `Travel.mat` memoizes `nloopcalc.m` and combines `TAV.mat` and `Codons.mat`, making it redundant. Huge performance increases resulted because `nloopcalc` was called so many times.
  * `dysentery.pl --rcheck`

### 3.14 ###
  * Back frameshifts are now detected. Every function that took a forward frameshift parameter now takes a back frameshift parameter as well.
  * `dysentery.pl` will convert amino acid primary structure sequences to codon sequences, randomly choosing codons. In addition, it verifies hand-customized codon sequences against a primary structure sequence.
  * Sine and cosine in calculating window function now bounded for more accurate results.
  * `jejunity.m` plots the average path taken by the codon in addition to error bars. It takes the same commands as `opportunity.m`, both displaying the same information about deviation.

### 3.1 ###
  * `opportunity.m`: Spits out superimposed displacement plots for rainy day family slideshows
  * `impunity.m`: Puts the yields for all sequences in a folder into a `results.txt` file
  * `unity.m` and `megaunity.m` and `cornerstone.pl` now web-enabled
  * `fabio.pl` (web-enabled): Parses TAV frequencies and spits them back out in the order of `Codons.mat`, i.e. spits out new `TAV.mat` for new species

### 3 ###
  * Unused functions deleted but available from [the archives](http://theframeshiftkids.googlecode.com/svn/tags/2+some/GSPtools/)
  * `fcalcmpx` split up and `megaunity` now several times faster
  * new `cornerstone.pl` parses Genbank search results and
  * no more fudging with free2bind (Thanks, Scott!)

## Why Version 1.1 ##
We re-released with Lalit's deterministic model merged other changes we've made since 3.1, merely replacing 3.1's `displacement.m`. Therefore it's web-enabled and has 3.1's Perl tools and performance.