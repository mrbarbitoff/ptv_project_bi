#!/usr/bin/env python3.6

from math import factorial
import pandas as pd
import numpy as np
#from scipy.stats import invgauss
from scipy.optimize import minimize
from scipy.integrate import quad

df = pd.read_csv('full_out.tsv', sep='\t')

# Get rid of zero in parameter fitting (change to None)
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
#    if res is None:
#        return 0
    return res

def P_n_integral(ig_params, df_sub):
    a, b = ig_params
    pars_zipped = zip(df_sub['AC_main'], df_sub['AN_main'], df_sub['U'])
    results = [quad(poisIG, 0, 1, args=(*x, a, b))[0] for x in pars_zipped]
    logL = np.array([np.log10(x) for x in results if x > 0])
    logL = logL.sum()
    print(f'alpha={a}, beta={b}, likelihood={logL}')
    return -logL

#def estimateS(n, N, U, a, b):
#    likelihoods = [(s, poisIG(s, n, N, U, a, b)) for s in np.linspace(0.0001, 1, 10000)]
#    S = likelihoods[int(np.argmax([x[1] for x in likelihoods if x[1] is not None]))]
#    return S

optimals = []
for tercile in range(1, 4):
    df_sub = df[df['tercile']== tercile]
    optres = minimize(P_n_integral, (1, 1), args=(df_sub),
        bounds=((0.001, 2), (0.0001, 2)), method="trust-constr") 
    optimals.append((optres.x, optres.fun))
    print(f'Optimal parameters {optres.x}, function value = {optres.fun} for tercile {tercile}')

print(optimals)

df['a'] = [optimals[x - 1][0][0] for x in df['tercile']]
df['b'] = [optimals[x - 1][0][1] for x in df['tercile']]

#input_params = zip(df['AC_main'], df['AN_main'], df['U'], df['a'], df['b'])

#S = []
#for gene in input_params:
#    S.append(estimateS(*gene))
#
#df['s_main'] = [x[0] for x in S]

df.to_csv('full_out_IG.tsv', index=False)
