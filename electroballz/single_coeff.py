from scipy import *
from scipy.special import sph_jn, sph_yn

# The following is an entirely computationally inefficient draft, intended for basic orientation.
# TODO clean up the code, possibly delete bf_coeff and mie_coeff_general, leave mie_coeff_mult only.

def jl(l,z):
    """Wrapper for sph_jn (discards the unnecessary data)"""
    return sph_jn(l, z)[0][l]

def jl_d(l,z):
    """Wrapper for sph_jn (discards the unnecessary data, keeps the first derivative)"""


def yl(l,z):
    """Wrapper for sph_yn (discards the unnecessary data)"""
    return sph_yn(l, z)[0][l]

# These are the Hankel equivalents of sph_jn, sph_yn
def sph_h1n(l,z):
    """Compute the spherical Hankel function h1(z) and its derivative for all orders up to and including n."""
    jn = sph_jn(l, z)
    yn = sph_yn(l, z)
    return (jn[0] + 1j*yn[0], jn[1] + 1j*yn[1])

def sph_h2n(l,z):
    """Compute the spherical Hankel function h2(z) and its derivative for all orders up to and including n."""
    jn = sph_jn(l, z)
    yn = sph_yn(l, z)
    return (jn[0] - 1j*yn[0], jn[1] - 1j*yn[1])

def h1l(l,z):
    """First spherical Hankel function"""
    return jl(l,z) + 1j*yl(l,z)

def h2l(l,z):
    """Second spherical Hankel function"""
    return j1(l,z) - 1j*yl(l,z)

# N.B. eta is sqrt(permeability/permittivity)


def bf_coeff(l, km, k0, etam, eta0, r): # THIS FUNCTION IS OBSOLETE.
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

# dictionary with available Bessel/Hankel functions
zfuncs = { 'j' : sph_jn, 'y' : sph_yn, 
    'h1' : sph_h1n, 'h2' : sph_h2n, 'h+' : sph_h1n, 'h-' : sph_h2n}

def mie_coeff_general(l, km, k0, etam, eta0, r, z_inc, z_out):
    """Scattering coefficients with optional external field base functions.
    
    Calculates Mie coefficients for arbitrary pair of "incoming" z_inc and
    "outcoming" z_out functions from the set of 'jl', 'yl', 'h1l', 'h2l'. 
    
    Returns a tuple with the ratios b/a and f/a, where
    a * v_inc is the "incoming" part of a spherical wave outside, of the given order l and which includes 
    the z_inc function; b * v_out is the corresponding "outcoming part". Inside the sphere, the solution
    is f * v, where v is the nonsingular (containing the Bessel function of first kind only) 
    spherical wave of the same order.

    """

    sph_z_inc = zfuncs[z_inc]
    sph_z_out = zfuncs[z_out]
    
    # The formulas are completely the same, so we just need to adjust the above code to the new names...
    sph_j_kmr = sph_jn(l, km*r)
    sph_z_inc_k0r = sph_z_inc(l, k0*r)
    sph_z_out_k0r = sph_z_out(l, k0*r)

    jm = sph_j_kmr[0][l]
    z0inc = sph_z_inc_k0r[0][l]
    z0out = sph_z_out_k0r[0][l]

    Jm = jm + km * r * sph_j_kmr[1][l]
    Z0inc = z0inc + k0 * r * sph_z_inc_k0r[1][l]
    Z0out = z0out + k0 * r * sph_z_out_k0r[1][l]

    denom1 = z0inc*Jm*k0*eta0 - Z0inc*jm*km*etam
    b1_a1 = - (z0out*Jm*k0*eta0 - Z0out*jm*km*etam) / denom1
    f1_a1 = - k0 * sqrt(eta0*etam) * (Z0inc*z0out - z0inc*Z0out) / denom1

    denom2 = Z0inc*jm*km*eta0 - z0inc*Jm*k0*etam
    b2_a2 = - (Z0out*jm*km*eta0 - z0out*Jm*k0*etam) / denom2
    f2_a2 = - k0 * sqrt(eta0*etam) * (-Z0inc*z0out + z0inc*Z0out) / denom2

    return (b1_a1, f1_a1, b2_a2, f2_a2)

def mie_coeff_mult(lmax, km, k0, etam, eta0, r, z_inc, z_out):
    """Scattering coefficients of order up to lmax with optional external field base functions.
    
    Calculates Mie coefficients for arbitrary pair of "incoming" z_inc and
    "outcoming" z_out functions from the set of 'jl', 'yl', 'h1l', 'h2l'. 
    
    Returns a tuple with the ratios b/a and f/a, where
    a * v_inc is the "incoming" part of a spherical wave outside, of the given order l and which includes 
    the z_inc function; b * v_out is the corresponding "outcoming part". Inside the sphere, the solution
    is f * v, where v is the nonsingular (containing the Bessel function of first kind only) 
    spherical wave of the same order.

    """

    sph_z_inc = zfuncs[z_inc]
    sph_z_out = zfuncs[z_out]
    
    # The formulas are completely the same, so we just need to adjust the above code to the new names...
    sph_j_kmr = sph_jn(lmax, km*r)
    sph_z_inc_k0r = sph_z_inc(lmax, k0*r)
    sph_z_out_k0r = sph_z_out(lmax, k0*r)

    jm = sph_j_kmr[0]
    z0inc = sph_z_inc_k0r[0]
    z0out = sph_z_out_k0r[0]

    Jm = jm + km * r * sph_j_kmr[1]
    Z0inc = z0inc + k0 * r * sph_z_inc_k0r[1]
    Z0out = z0out + k0 * r * sph_z_out_k0r[1]

    denom1 = z0inc*Jm*k0*eta0 - Z0inc*jm*km*etam
    b1_a1 = - (z0out*Jm*k0*eta0 - Z0out*jm*km*etam) / denom1
    f1_a1 = - k0 * sqrt(eta0*etam) * (Z0inc*z0out - z0inc*Z0out) / denom1

    denom2 = Z0inc*jm*km*eta0 - z0inc*Jm*k0*etam
    b2_a2 = - (Z0out*jm*km*eta0 - z0out*Jm*k0*etam) / denom2
    f2_a2 = - k0 * sqrt(eta0*etam) * (-Z0inc*z0out + z0inc*Z0out) / denom2

    return (b1_a1, f1_a1, b2_a2, f2_a2)


