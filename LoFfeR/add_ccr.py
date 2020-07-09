#!/usr/bin/env python3.6

import sys
import numpy as np

varfile, ccrfile = sys.argv[1:]

var_ccrs = {}
with open(ccrfile, 'r') as ccrf:
    for line in ccrf:
        content = line.strip().split('\t')
        pos = int(content[1]) + 3
        var_id = f'{content[0]}_{pos}'
        var_ccr = var_ccrs.get(var_id, [0, 0, 0])
        if float(content[7]) > var_ccr[0]:
            var_ccrs[var_id] = [float(x) for x in content[7:10]]

#with open(varfile, 'r') as vfile:
#    for line in vfile:
#        print(f'{line.strip()}\tccr_exon_pct\tccr_exon_z\tccr_exon_frac')
#        break

with open(varfile, 'r') as vfile:
    for line in vfile:
        if line.startswith('Chr') or 'REF\tALT' in line:
            print(f'{line.strip()}\tccr-exon_pct\tccr_exon_z\tccr_exon_frac')
            continue
        content = line.strip().split('\t')
        var_id = '_'.join(content[:2])
        if var_id in var_ccrs:
            to_out = content + [str(x) for x in var_ccrs[var_id]]
        else:
            to_out = content + ['NA', 'NA', 'NA']
        print('\t'.join(to_out))
