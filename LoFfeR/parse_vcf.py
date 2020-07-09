#!/usr/bin/env python3.6

import numpy as np
import pandas as pd
import gzip
import re
import sys

vcf_path, gtex_file = sys.argv[1:]

with open(vcf_path, mode='rt') as f:
        for line in f:
                line = line.rstrip().split('\t')
                if '#' not in line[0]:
                        ALT=line[4]
                        INFO=line[7].split(';')
                        break

# pattern = '^AN$|^AC$|^AC_\w{3}$|^AN_\w{3}$|^vep'
pattern = '^CSQ'
fields = [key for key in [i.split('=')[0] for i in INFO] if re.search(pattern, key) != None]
# fields = [i for i in fields if all(['oth' not in i,
#                                                                       'raw' not in i])]

# info_keys = []
# for k in fields:
#       if k != 'vep':
#               info_keys.append(k)
                
df_tr = pd.read_csv(gtex_file, '\t')
df_tr['gene_id'] = df_tr['gene_id'].str.split('.').str[0]
df_tr['transcript_id'] = df_tr['transcript_id'].str.split('.').str[0]
df_tr.head()

def process_line(line, fields): #, info_keys):
        
        line = line.split('\t')
        INFO=line[7].split(';')
        if 'CSQ=' not in line[7] or len(line) > 10:
            return None
#       pattern = '^AN$|^AC$|^AC_\w{3}$|^AN_\w{3}$|^vep'
        pattern = '^CSQ'
        fields = [key for key in [i.split('=')[0] for i in INFO] if re.search(pattern, key) != None]
#       fields = [i for i in fields if all(['oth' not in i,
#                                                                               'raw' not in i])]
        
        INFO_index = sum([[i for i, item in enumerate(INFO) if re.search(f'^{regex}=', item)] for regex in fields], [])
        INFO_index = dict(zip(fields, INFO_index))
#        print(INFO_index)

        tr_list = (INFO[INFO_index['CSQ']].split('=')[1]).split(',')
        
        for i in tr_list:

                col_dict = {}
        
                col_dict['Chr'] = []
                col_dict['pos'] = []
                col_dict['rs_ID'] = []
                
                col_dict['is_mean'] = []
                col_dict['is_max'] = []
                col_dict['is_5pct_mean'] = []
                col_dict['is_5pct_max'] = []
                col_dict['median_mean'] = []
                
#               col_dict = {**col_dict,
#                                       **{k:v for k, v in dict(zip(fields, [[] for i in range(len(fields))])).items() if 'vep' not in k}}
           
                col_dict['REF'] = []
                col_dict['ALT'] = []
                col_dict['gene'] = []
                col_dict['conseq'] = []
                col_dict['ens_tr'] = []

                spl = i.split('|')
                
                conseq, gene_name, canonical, LoF_filter, LoF, ens_tr = spl[1], spl[3], spl[23], spl[-3], spl[-2], spl[6]
#                print(f'{conseq} {gene_name} {canonical} {LoF_filter} {LoF} {ens_tr}')

                article_filter_success = False

                if all([any(["stop_gained" in conseq, "splice_donor" in conseq, "splice_acceptor" in conseq]), LoF_filter == '', LoF == '']):
#                        print('Got here')
                        for k, v in col_dict.items():
                        
                                if all(['REF' not in k,
                                                'ALT' not in k,
                                                'gene' not in k,
                                                'ens_tr' not in k,
                                                'Chr' not in k,
                                                'pos' not in k,
                                                'rs_ID' not in k,
                                                'conseq' not in k,
                                                'is_mean' not in k,
                                                'is_max' not in k,
                                                'is_5pct_mean' not in k,
                                                'is_5pct_max' not in k,
                                                'median_mean' not in k]):
                                
                                        v.append(int(INFO[INFO_index[k]].split('=')[1]))
                        
                        col_dict['Chr'].append(line[0])       
                        col_dict['pos'].append(line[1])
                        col_dict['rs_ID'].append(line[2])
                        col_dict['gene'].append(gene_name)
                        col_dict['ens_tr'].append(ens_tr)
                        col_dict['REF'].append(line[3])
                        col_dict['ALT'].append(line[4])
                        col_dict['conseq'].append(conseq)
                        
                        article_filter_success = True

                else:
                        continue

                if article_filter_success:
                                
                        for i in tr_list:
                                
                                conseq, tr = i.split('|')[1], i.split('|')[6]
                                df = df_tr[df_tr['transcript_id']==tr]
                                
                                if df.shape[0]==0:
                                        continue
                                        
                                if any(["stop_gained" in conseq, "splice_donor" in conseq, "splice_acceptor" in conseq]):

                                                col_dict['is_mean'].append(df['is_max_mean'].iloc[0])
                                                col_dict['is_max'].append(df['is_max_max'].iloc[0])
                                                col_dict['is_5pct_mean'].append(df['is_5pct_max_mean'].iloc[0])
                                                col_dict['is_5pct_max'].append(df['is_5pct_max_max'].iloc[0])
                                                col_dict['median_mean'].append(float(df['max_mean'].iloc[0]))

                        if len(col_dict['is_mean'])>0:

                                col_dict['is_mean'] = [1] if 1 in col_dict['is_mean'] else [0]
                                col_dict['is_max'] = [1] if 1 in col_dict['is_max'] else [0]
                                col_dict['is_5pct_mean'] = [1] if 1 in col_dict['is_5pct_mean'] else [0]
                                col_dict['is_5pct_max'] = [1] if 1 in col_dict['is_5pct_max'] else [0]
                                col_dict['median_mean'] = [np.nanmedian(col_dict['median_mean'])]
                                
                        else:
                                col_dict['is_mean'] = [np.nan]
                                col_dict['is_max'] = [np.nan]
                                col_dict['is_5pct_mean'] = [np.nan]
                                col_dict['is_5pct_max'] = [np.nan]
                                col_dict['median_mean'] = [np.nan]
                        
                        if len(col_dict['gene']) == 1:
                                
                                return col_dict
                                
                        else:
                                col_dict['Chr'] = col_dict['Chr'][0]
                                col_dict['pos'] = col_dict['pos'][0]
                                col_dict['rs_ID'] = col_dict['rs_ID'][0]
                                col_dict['REF'] = col_dict['REF'][0]
                                col_dict['ALT'] = col_dict['ALT'][0]
                                
                                for k, v in col_dict.items():
                                        
                                        if all(['ALT' not in k,
                                                        'gene' not in k,
                                                        'ens_tr' not in k,
                                                        'Chr' not in k,
                                                        'pos' not in k,
                                                        'rs_ID' not in k,
                                                        'conseq' not in k,
                                                        'is_mean' not in k,
                                                        'is_max' not in k,
                                                        'is_5pct_mean' not in k,
                                                        'is_5pct_max' not in k,
                                                        'median_mean' not in k]):
                                                
                                                col_dict[k] = v[0]
                                                
                                col_dict['gene'] = ','.join(col_dict['gene'])
                                col_dict['conseq'] = ','.join(col_dict['conseq'])
                                col_dict['ens_tr'] = ','.join(col_dict['ens_tr'])
                                
                                return col_dict

with open('ptv_data.tsv', 'w') as w:
        w.write('Chr\tpos\trs_ID\tis_mean\tis_max\tis_5pct_mean\tis_5pct_max\tmedian_mean\tREF\tALT\tgene\tconseq\tens_tr\n')

with open(vcf_path, 'rt') as f:
        for line in f:
                if '#' not in line:
                        if line.split()[6] in ['PASS', '.', '']:

                                d = process_line(line, fields)
                                if d:
                                        pass
#                                        print('Emitted')
                                        d = pd.DataFrame(d)
                                        d.to_csv('ptv_data.tsv', sep='\t', mode='a', header=None, index=None)


