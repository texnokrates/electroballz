# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from electroballz.translations import translation_S_real
from electroballz.single_coeff import jl, h1l
import scipy.constants as sc

# <codecell>

N=10
k=2*sc.pi/(500*sc.nano)
d=100*sc.nano
lp=1
mp=1
sigmap=0
taup=1
eta = 0 # Zkusit i pi/2
psi = 0
z = jl

for l in range(0,N):
    for m in range (0,l+1):
        for sigma in ['e','o']:
            for tau in [1,2]:
                print(l,m,sigma,tau,
                      translation_S_real(tau, sigma, m, l, d, taup, sigmap, mp, lp, eta, psi, z, k))
                
