#!/usr/bin/env python3.6

import gzip
import sys
import re

vars, pext = sys.argv[1:]

pexts = {}
with gzip.open(pext, 'rt') as pext_file:
    for line in pext_file:
        if line.startswith('ensg'):
            continue
        content = line.strip().split('\t')
        if content[3] == 'NaN':
            continue
        pexts[content[1]] = content[3]

with open(vars, 'r') as vfile:
    for line in vfile:
        if line.startswith('Chr'):
            print(f'{line.strip()}\tpext_avg')
            continue
        content = line.strip().split('\t')
        pos = ':'.join(content[:2])
        pext_avg = 'NaN'
        if pos not in pexts:
            for i in range(-3, 4):
                pos = content[0] + ':' + str(int(content[1]) + i)
                if pos not in pexts:
                    continue
                else:
                    pext_avg = pexts[pos]
                    break
        else:
            pext_avg = pexts[pos]
        print(f'{line.strip()}\t{pext_avg}')
