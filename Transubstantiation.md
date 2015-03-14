A note on lexicon: ++r means with `reloop` and ++s means with `sofar/acc`. The minus signs mean without.

Without loss of generality, we'll just consider going from ++r++s to --r--s. All these changes go into `displacement.m`.

  * Add a line to set `acc` in `probabilities` to "1" on the very first line.
  * Change `age_limit` to 1000.
  * Change the `1` in the `reloop` line to `0`.