#!/usr/bin/env python3.6

import re
import sys

vcf, vardata, constraint, variant_id_cols, out_file = sys.argv[1:]
variant_id_cols = [int(x) for x in variant_id_cols.split(',')]

constraint_measures = {}
with open(constraint, 'r') as cf:
    for line in cf:
        if line.startswith('gene'):
            continue
        content = line.strip().split('\t')
        transcript = content[1]
#        print(transcript, content[21], content[30])
        constraint_measures[transcript] = (content[5], content[21], content[30])

variant_stats = {}
with open(vcf, 'r') as vfile:
    for line in vfile:
        if line.startswith('#'):
            continue
        if 'CSQ=' in line.split('\t')[7]:
            splitter_field = 'CSQ='
        elif 'vep=' in line.split('\t')[7]:
            splitter_field = 'vep='
        else:
            continue
        content = line.strip().split('\t')
        variant_id = ':'.join(content[:5])
#        print(variant_id)
#        if variant_id in variant_stats:
#            print('Strange duplicate')
        info_field = content[7]
#        print(info_field)
        effects = info_field.split(splitter_field)[1].split(';')[0].split(',')
        for effect in effects:
            annotations = effect.split('|')
            if annotations[0] == content[4]:
                tr_id = annotations[6]
                tr_stats = constraint_measures.get(tr_id, ['NA', 'NA', 'NA'])
                if annotations[1] == 'stop_gained' or 'splice_' in annotations[1]:
                    var_stats = variant_stats.get(variant_id, [0, 1, False, 1, 1])
                    if any([tr_stats[1] == 'NA', tr_stats[2] == 'NA']):
                        variant_stats[variant_id] = var_stats
                        continue
                    if float(tr_stats[1]) > var_stats[0]:
                        var_stats[0] = float(tr_stats[1])
                    if float(tr_stats[2]) < var_stats[1]:
                        var_stats[1] = float(tr_stats[2])
                    variant_stats[variant_id] = var_stats
                else:
                    var_stats = variant_stats.get(variant_id, [0, 1, False, 1, 1])
                    var_stats[2] = True
                    if any([tr_stats[1] == 'NA', tr_stats[2] == 'NA']):
                        variant_stats[variant_id] = var_stats
                        continue
                    if float(tr_stats[0]) < var_stats[3]:
                        var_stats[3] = float(tr_stats[0])
                    if float(tr_stats[2]) < var_stats[4]:
                        var_stats[4] = float(tr_stats[2])
                    variant_stats[variant_id] = var_stats
        if variant_stats.get(variant_id, 'xyi') == 'xyi':
            variant_stats[variant_id] = ['NA', 'NA', 'NA', 'NA', 'NA']


#print(len(variant_stats))

counter = 0
with open(vardata, 'r') as vdf, open(out_file, 'w') as plif:
    for line in vdf:
        if line.startswith('Chr') or 'REF\tALT' in line:
            plif.write(line.strip() + '\tmax_pLI\tmin_LOEUF\trescue\trescue_oe_mis\trescue_loeuf_min\n')
            continue
        content = line.strip().split('\t')
        var_id = ':'.join([content[x] for x in variant_id_cols])
        if var_id in variant_stats and variant_stats[var_id][0] != 'NA':
            counter += 1
 #       print(var_id)
        var_stats = [str(x) for x in variant_stats.get(var_id, ['NA', 'NA', 'NA', 'NA', 'NA'])]
        plif.write(line.strip() + '\t' + '\t'.join(var_stats) + '\n')

#print(counter)

