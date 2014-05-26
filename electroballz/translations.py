import clebsch from qutip

#TODO check data types, optimize
#!!!! Je clebsch ta správná funkce pro 3j symbol?

def wigner3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol"""
    # see http://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients#Relation_to_3-jm_symbols
    return clebsch(j1, j2, j3, m1, m2, -m3) * (-1)**(j1-j2-m3) / (
            sqrt(2*j3 + 1))




def factorial_frac(numer, denom):
    """Calculates numer!/denom!"""
    # we silently suppose that the arguments are integer
    # slow&dirty way for now (fix that later)
    return numer / float(denom)


def bks_C(m, l, mp, lp, d, eta, z):
    """[BKS 1991], (5.36)"""
    coseta = cos(eta)
    res = 0
    for lambd in range(abs(l-lp) + 1, l + lp + 1):
        # If l + l' + λ is od, the first 3j symbol is zero.
        if 0 == (l + lp + lambd) % 2:
            line2 = (-1)**(lp - l + lambd) * (2*lambd + 1) * sqrt(
                    #TODO podíl faktorů udělat nějak chytřeji. Co dělá python s přetečením?
                    (2*l+1)*(2*lp+1)*factorial_frac(lambd-m+mp,lambd+m-mp) /
                    (l * (l+1) * lp * (lp + 1)))
            line3 = wigner3j(l, lp, lambd, 0, 0, 0)
            line3 *= wigner3j(l, lp, lambd, m, -mp, mp - m)
            line3 *= l*(l+1) + lp*(lp+1) - lambd*(lambd+1)


             
            

    

