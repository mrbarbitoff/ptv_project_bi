#!/usr/bin/env python3.6

from math import factorial
import pandas as pd
import numpy as np
#from scipy.stats import invgauss
#from scipy.optimize import minimize
#from scipy.integrate import quad
from scipy.stats import poisson
from concurrent.futures import ProcessPoolExecutor

df_full = pd.read_csv('full_corr_wshet.tsv')

def poisIG(s, n, N, U, a, b):
    nu = N * U
    lamb = nu/s
    try:
        pois_enum = lamb**n * np.exp(-lamb)
        pois = pois_enum/factorial(n)
    except OverflowError:
        normp_fc = 1/((2 * np.pi * lamb)**0.5)
        normp_sc = np.exp(-(((n - lamb)**2)/(2 * lamb)))
        pois = normp_fc * normp_sc
    fc_ig = (b/(2 * np.pi * (s**3)))**0.5
    sc_ig = np.exp(-((b * ((s - a)**2))/(2 * (a**2) * s)))
    ig = fc_ig * sc_ig
    res = pois*ig if abs(pois * ig) != np.inf else None
    if res is None:
        return 0
    return res

def P_n(n, N, U, a, b):
#    print(n, N, U, a, b)
    result = np.sum([poisIG(s, n, N, U, a, b) * 0.0001 for s in np.linspace(0.0001, 1, 10000)])
    return result

def calculate_CPL(input_values):
    pop_n, pop_N, U, a, b = input_values
    denoms = np.array([P_n(*vals, U, a, b) for vals in zip(pop_n, pop_N)])
#    print(denoms)
    def enums(s, pop_n, pop_N, U, a, b):
        return np.array([poisIG(s, *vals, U, a, b) * 0.0001 for vals in zip(pop_n, pop_N)])
    densities = []
    for s in np.linspace(0.0001, 1, 10000):
        enum = enums(s, pop_n, pop_N, U, a, b)
        parts = np.array(enum/denoms)
        if np.sum(enum > denoms):
            print('ERROR')
        densities.append(np.prod(parts))
    if np.sum(densities) > 1:
        print((pop_n, pop_N, U, a, b, denoms, densities))
    return -np.log10(np.sum(densities))

def calculate_replicate(input_list):
    return list([calculate_CPL(input_values) for input_values in input_list])

ns = list(zip(df_full['AC_afr'], df_full['AC_amr'], df_full['AC_eas'], df_full['AC_nfe'], df_full['AC_sas']))
Ns = list(zip(df_full['AN_afr'], df_full['AN_amr'], df_full['AN_eas'], df_full['AN_nfe'], df_full['AN_sas']))
in_values = list(zip(ns, Ns, df_full['U'], df_full['a'], df_full['b']))
print('Starting true CPL calculation...')

CPLs = [calculate_CPL(inputs) for inputs in in_values]
df_full['L'] = CPLs
print('True values calculated, proceeding to replicates...')

nus = np.array([list(y * x[1] for y in x[0]) for x in zip(Ns, df_full['U'])])
mus = list([list(y/x[1] for y in x[0]) for x in zip(nus, df_full['s_main'])])
exps = list([list([list([poisson.rvs(x) for x in y]) for y in mus]) for i in range(100)])
inputs = list([list(zip(x, Ns, df_full['U'], df_full['a'], df_full['b'])) for x in exps])
print('Input values generated succesfully, proceeding...')

with ProcessPoolExecutor(max_workers=34) as threads:
    replicates = threads.map(calculate_replicate, inputs)

print('Done parallel execution, calculating summaries and writing...')

replicates_df = pd.DataFrame(replicates)
df_full['L_mean'] = list(replicates_df.apply(np.mean, axis=0))
df_full['L_sd'] = list(replicates_df.apply(np.std, axis=0))
df_full['L_replicate'] = replicates_df.loc[0]

df_full.to_csv('full_out_boot_stats.csv')
