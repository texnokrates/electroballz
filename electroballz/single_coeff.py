from scipy import *
from scipy.special import sph_jn, sph_yn

# The following is an entirely computationally inefficient draft, intended for basic orientation.

def jl(l,z):
    """Wrapper for sph_jn (discards the unnecessary data)"""
    return sph_jn(l, z)[0][l]

def yl(l,z):
    """Wrapper for sph_yn (discards the unnecessary data)"""
    return sph_yn(l, z)[0][l]

def h1l(l,z):
    """First spherical Hankel function"""
    return jl(l,z) + 1j*yl(l,z)

def h2l(l,z):
    """Second spherical Hankel function"""
    return j1(l,z) - 1j*yl(l,z)

def bf_coeff(l, km, k0, etam, eta0, r):
    """Ratios between (b1lm,f1lm) and a1lm. See the single_spherical_wave_scatter.nb file"""
    sph_j_kmr = sph_jn(l, km*r)
    sph_j_k0r = sph_jn(l, k0*r)
    sph_y_k0r = sph_yn(l, k0*r)

    jm = sph_j_kmr[0][l]
    h01 = sph_j_k0r[0][l] + 1j * sph_y_k0r[0][l]
    h02 = sph_j_k0r[0][l] - 1j * sph_y_k0r[0][l]

    Jm = jm + km * r * sph_j_kmr[1][l]
    H01 = h01 + k0 * r * (sph_j_k0r[1][l] + 1j * sph_y_k0r[1][l])
    H02 = h02 + k0 * r * (sph_j_k0r[1][l] - 1j * sph_y_k0r[1][l])

    denom1 = h01*Jm*k0*eta0 - H01*jm*km*etam
    b1_a1 = - (h02*Jm*k0*eta0 - H02*jm*km*etam) / denom1
    f1_a1 = - k0 * sqrt(eta0*etam) * (H01*h02 - h01*H02) / denom1

    denom2 = (H01*jm*km*eta0 - h01*Jm*k0*etam)
    b2_a2 = - (H02*jm*km*eta0 - h02*Jm*k0*etam) / denom2
    f2_a2 = - k0 * sqrt(eta0*etam) * (-H01*h02 + h01*H02) / denom2
  
    return (b1_a1, f1_a1, b2_a2, f2_a2)


