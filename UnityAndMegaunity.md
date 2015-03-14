## Introduction ##
Currently, this page is just a series of emails that sheds light on what unity/megaunity do. (Editing it into actual documentation is welcome.)

## Email 1 ##
> unity graphs the displacement of the reading frame from the expected
> reading frame (so a plot that goes from 0 to 2 indicates a ~~+2~~ +1
> frameshift). megaunity simply repeats unity, tabulating the number of
> times it frameshifts at codon uga-25 in prfB without
> mistakes, thus calculating a probability and therefore a "yield" of
> sorts.

> megaunity also outputs all the frameshifts it detects; however Dr.
> Bitzer has made it clear that only the uga-25 frameshift is correct.
> Any other frameshift is a mistake made by the ribosome.

> megaunity isn't useful unless you're running prfB
> through it or you know the correct frameshifts ahead of time and can
> modify the program to check for those.

## Email 2 ##

> Just to clarify what Hao wrote a bit earlier, unity will input a genetic sequence (with a 12-base leader and ending at the stop codon) and output a displacement plot.  A jump from 0 to 2 indicates a +1 frameshift.  (Hao wrote a +2 frameshift by accident, but in actuality we have not changed the method of drawing the graph, except for incoporating probability.)

> The uga-25 frameshift is the only one that is actually part of prfB, obviously.  All other frameshifts are unwanted, when the ribsome just chooses the wrong codon.  Please correct us if we are wrong, but we assume that if the ribosome chooses the wrong codon, it continues in the wrong frame (this is what we think Lalit's program used to do, so we kept that feature in this program).

> Megaunity does what Hao wrote.  If you are testing a gene that does not frameshift, you can just go in and modify the line that compares a variable (called "ants", I think), to "uga; 25"; change it to compare to a blank string.
