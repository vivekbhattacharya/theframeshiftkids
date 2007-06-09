lines = open('tav.txt', 'r').readlines()
a = []

import re
splitter = re.compile(r'\s+')
for line in lines:
    line = line.strip()
    if not line: continue
    if line.startswith('C'): continue
    a += splitter.split(line)

p = open('tav2.txt', 'w+')
p.writelines('\n'.join(a))