## Introduction ##

We propose that "bumps" in a gene's displacement plot affect frameshifts later in the graph. Thus, "smoothing out" a graph from the start will reduce the number of frameshifts predicted by the model, i.e. the number of wrong codons chosen.


## Outline ##
  * Identifies areas with large amounts of frameshifting, i.e. volatile areas where the ribosome does not choose the right codon often.
  * Selects an area four codons 5' of the target codon.
  * Randomly permutes with synonyms of these four codons to minimize frameshifting at or near that given codon. However, it does not try to optimize _yield_ by permuting codons.
  * Repeats the process for multiple volatile areas.

Thrushbaby is especially effective if run multiple times, each time on the improved sequence.