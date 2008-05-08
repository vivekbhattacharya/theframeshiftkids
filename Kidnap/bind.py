from freier import Freier
from strand import Strand

class BindSequenceException(Exception): pass
class Bind(object):
    def __init__(self, cat = Freier(273.15+37)):
        self.best_score = 0
        self.doublet = cat

    def energy(self, seq, rna):
        length = len(rna)
        if not length == len(seq):
            raise BindSequenceException('Unequal sequence lengths')
        if length < 2:
            raise BindSequenceException('Sequences are too short')

        best_score = score = helix_end = 0

        # I need at least two bases in each sequence in order to form
        # a structure between them.
        x = Strand(seq, length)
        y = Strand(rna, length)

        # Optimize away dot access
        di, dt = self.doublet.internal, self.doublet.terminal

        for i in xrange(0, length - 1):
            x.update(i); y.update(i)

            # Start each helix structure with a terminal doublet.
            if score == 0: score = dt(x.all, y.all, True)
            else:
                # Some internal "doublets" need more context to
                # be scored correctly.  See the scoring for GU
                # in XiaMathews#internal_doublet for more details.
                score += di(x.all, y.all, x.context, y.context)

            if score > 0: score = 0
            elif score < best_score:
                best_score = score
                helix_end = i

        # Swap out the last internal for a terminal doublet score
        # so that the helix starts and ends with a terminal doublet.
        # We have to do this now because we cannot know in advance
        # where this terminal doublet is going to be.
        #   We should check for "self symmetry." With the force_bind,
        # we are assuming that the strands are stretched out, and thus
        # not stuck to themselves forming their own hairpin loop. The
        # "self symmetry" penalty arises from the individual strands
        # being stuck together in their own hairpins.
        if best_score < 0:
            x.update(helix_end); y.update(helix_end)

            internal_score = di(x.all, y.all)
            terminal_score = dt(x.all, y.all, False)
            best_score -= internal_score - terminal_score

        return best_score + self.doublet.init_penalty
