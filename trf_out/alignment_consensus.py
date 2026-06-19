#!/usr/bin/env python3
import sys
from collections import Counter

seqs = []
name = None
buf = []

for line in open(sys.argv[1]):
    line = line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if buf:
            seqs.append("".join(buf).upper())
        buf = []
    else:
        buf.append(line)

if buf:
    seqs.append("".join(buf).upper())

L = max(map(len, seqs))
seqs = [s.ljust(L, "-") for s in seqs]

cons = []
for i in range(L):
    col = [s[i] for s in seqs if s[i] not in "-N?"]
    if not col:
        continue
    base, count = Counter(col).most_common(1)[0]
    cons.append(base)

print(">consensus")
print("".join(cons))
