class Strand(object):
    def __init__(self, seq, length):
        self.seq = list('*' + seq + '*')
        self.length = length

    def update(self, i):
        (s55, s5, s3, s33) = self.seq[i:i+4]
        self.all = (s5, s3)
        
        self.context = ['', '']
        if i > 0 and i < self.length - 2:
            self.context = (s55, s33)