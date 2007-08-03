""" Turning values from Freier 1986 to code """

import sys
sys.path = ['..'] + sys.path

BIGNUM = 424242
from sequence import valid_pair, wc_pair

class Freier(object):
    def __init__(self, temp):
        self.temp = temp
        self.init_penalty = 3.4

    # Get rid of context
    def internal(self, top, bottom, *args):
        t5, t3, b3, b5 = top + bottom
        if not(valid_pair(t5, b3) and valid_pair(t3, b5)): return BIGNUM
        
        scores = {
            'au': {
                # Watson-Crick matches and G/U mismatches
                'au': -0.9, 'ua': -0.9, 'gc': -1.7, 'cg': -2.1,
                'gu': -0.5, 'ug': -0.7,
            }, 'ua': {
                'au': -1.1, 'ua': -0.9, 'gc': -1.8, 'cg': -2.3,
                'gu': -0.7, 'ug': -0.5,
            }, 'cg': {
                'au': -1.8, 'ua': -1.7, 'gc': -2.9, 'cg': -2.0,
                'gu': -1.5, 'ug': -1.5,
            }, 'gc': {
                'au': -2.3, 'ua': -2.1, 'gc': -2.9, 'cg': -3.4,
                'gu': -1.3, 'ug': -1.9,
            }, 'gu': {
                'au': -0.5, 'ua': -0.7, 'gc': -1.5, 'cg': -1.9,
                'gu': -0.5, 'ug': -0.5,
            }, 'ug': {
                'au': -0.7, 'ua': -0.5, 'gc': -1.5, 'cg': -1.3,
                'gu': -0.6, 'ug': -0.5,
            },
        }
        
        return scores[t5 + b3][t3 + b5]
    
    def terminal(self, top, bottom, left_side):
        t5, t3, b3, b5 = top + bottom
        if left_side: (t5, b5) = (b5, t5); (t3, b3) = (b3, t3)
        
        if not wc_pair(t5, b3): return BIGNUM
        scores = {
            'au': {
                # Watson-Crick matches and G/U mismatches
                'au': -0.9, 'ua': -0.9, 'gc': -1.7, 'cg': -2.1,
                'gu': -0.9, 'ug': -0.9,
                # Mismatches
                'aa': -0.8, 'cc': -0.7, 'gg': -1.0, 'uu': -0.8,
                'ac': -1.0, 'ag': -1.0, 'ca': -0.7, 'ga': -0.8,
                'uc': -0.8, 'cu': -0.7,
            }, 'ua': {
                'au': -1.1, 'ua': -0.9, 'gc': -1.8, 'cg': -2.3,
                'gu': -0.9, 'ug': -1.0,
                'aa': -1.0, 'cc': -0.6, 'gg': -1.2, 'uu': -0.5,
                'ac': -0.8, 'ag': -1.1, 'ca': -0.7, 'ga': -1.1,
                'uc': -0.6, 'cu': -0.5,
            }, 'cg': {
                'au': -1.8, 'ua': -1.7, 'gc': -2.0, 'cg': -2.9,
                'gu': -1.6, 'ug': -1.9,
                'aa': -1.9, 'cc': -1.1, 'gg': -1.9, 'uu': -1.2,
                'ac': -2.0, 'ag': -1.9, 'ca': -1.0, 'ga': -1.0,
                'uc': -1.5, 'cu': -0.8,
            }, 'gc': {
                'au': -2.3, 'ua': -2.1, 'gc': -2.9, 'cg': -3.4,
                'gu': -1.4, 'ug': -2.3,
                'aa': -1.1, 'cc': -0.6, 'gg': -1.4, 'uu': -0.7,
                'ac': -1.7, 'ag': -1.3, 'ca': -1.1, 'ga': -1.6,
                'uc': -0.8, 'cu': -0.5,
            },
        }
        
        return scores[t5 + b3][t3 + b5]