"""
Dielectric properties of wet soil.

Hallikainen et al. (1985) empirical model for complex relative dielectric
constant as a function of volumetric soil moisture, frequency, and texture.

De Zan et al. (2014) Eq. (5): vertical complex wavenumber in soil.
"""
import numpy as np

# Hallikainen 1985 coefficients, visually checked against MIMICS Table E.1,
# https://web.eecs.umich.edu/~lep/mimics.pdf, printed E-3/E-4 (PDF 146/147).
# Corrected 2026-10-15: the former 5 GHz row was mislabeled 4 GHz and its
# imaginary quadratic row was incorrect. Historical simulations must be rerun.
# eps_r = (a0 + a1*S + a2*C) + (b0 + b1*S + b2*C)*mv + (c0 + c1*S + c2*C)*mv^2
# S,C in the polynomial are PERCENT BY WEIGHT; API inputs remain fractions.
_HALLIKAINEN_COEFFS = {
    # freq_GHz: (real_coeffs[3x3], imag_coeffs[3x3])
    1.4: {
        "real": np.array([
            [2.862, -0.012,  0.001],
            [3.803,  0.462, -0.341],
            [119.006, -0.500,  0.633],
        ]),
        "imag": np.array([
            [0.356, -0.003, -0.008],
            [5.507,  0.044, -0.002],
            [17.753, -0.313,  0.206],
        ]),
    },
    4.0: {
        "real": np.array([
            [2.927, -0.012, -0.001],
            [5.505,  0.371,  0.062],
            [114.826, -0.389, -0.547],
        ]),
        "imag": np.array([
            [0.004,  0.001,  0.002],
            [0.951,  0.005, -0.010],
            [16.759,  0.192,  0.290],
        ]),
    },
}

for _freq, _real, _imag in [
    (6., [1.993,.002,.015,38.086,-.176,-.633,10.720,1.256,1.522],
         [-.123,.002,.003,7.502,-.058,-.116,2.942,.452,.543]),
    (8., [1.997,.002,.018,25.579,-.017,-.412,39.793,.723,.941],
         [-.201,.003,.003,11.266,-.085,-.155,.194,.584,.581]),
    (10., [2.502,-.003,-.003,10.101,.221,-.004,77.482,-.061,-.135],
          [-.070,.000,.001,6.620,.015,-.081,21.578,.293,.332]),
    (12., [2.200,-.001,.012,26.473,.013,-.523,34.333,.284,1.062],
          [-.142,.001,.003,11.868,-.059,-.225,7.817,.570,.801]),
    (14., [2.301,.001,.009,17.918,.084,-.282,50.149,.012,.387],
          [-.096,.001,.002,8.583,-.005,-.153,28.707,.297,.357]),
    (16., [2.237,.002,.009,15.505,.076,-.217,48.260,.168,.289],
          [-.027,-.001,.003,6.179,.074,-.086,34.126,.143,.206]),
    (18., [1.912,.007,.021,29.123,-.190,-.545,6.960,.822,1.195],
          [-.071,.000,.003,6.938,.029,-.128,29.945,.275,.377]),
]:
    _HALLIKAINEN_COEFFS[_freq] = dict(real=np.array(_real).reshape(3,3),
                                    imag=np.array(_imag).reshape(3,3))


def hallikainen_dielectric(mv, freq_ghz=1.4, sand=0.51, clay=0.13):
    """
    Complex relative dielectric constant of wet soil.

    Hallikainen et al. (1985), IEEE TGRS, Table 1.
    Reproduced as Fig. 2 in De Zan et al. (2014).

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture [0, 0.5].
    freq_ghz : float
        Radar frequency in GHz, 1.4–18. Linear interpolation between adjacent
        tabulated frequencies; no extrapolation. 5.405 uses the 4 and 6 GHz rows.
    sand : float
        Sand fraction [0, 1].
    clay : float
        Clay fraction [0, 1].

    Returns
    -------
    eps : complex ndarray
        Complex relative dielectric constant.
    """
    mv = np.asarray(mv, dtype=float)
    sand, clay = np.broadcast_arrays(np.asarray(sand, float), np.asarray(clay, float))
    if (not np.isfinite(freq_ghz) or not 1.4 <= freq_ghz <= 18
            or not np.all(np.isfinite(sand)) or not np.all(np.isfinite(clay))
            or np.any(sand < 0) or np.any(clay < 0) or np.any(sand+clay > 1)):
        raise ValueError('Frequency must be 1.4–18 GHz and texture fractions sum to <=1')
    if freq_ghz not in _HALLIKAINEN_COEFFS:
        # linear interpolation between nearest supported frequencies
        freqs = sorted(_HALLIKAINEN_COEFFS.keys())
        j = np.searchsorted(freqs, freq_ghz)
        f0, f1 = freqs[j-1], freqs[j]
        t = (freq_ghz - f0) / (f1 - f0)
        eps0 = hallikainen_dielectric(mv, f0, sand, clay)
        eps1 = hallikainen_dielectric(mv, f1, sand, clay)
        return eps0 + t * (eps1 - eps0)

    c = _HALLIKAINEN_COEFFS[freq_ghz]
    # Broadcast texture maps with moisture.  The previous implementation
    # created a 1-D array and therefore failed for per-pixel 2-D sand/clay.
    S, C = 100*sand, 100*clay
    texture = np.stack((np.ones_like(S), S, C), axis=-1)

    def _poly(coeffs):
        a = np.einsum("...j,ij->...i", texture, coeffs)
        return a[..., 0] + a[..., 1] * mv + a[..., 2] * mv**2

    eps_real = _poly(c["real"])
    eps_imag = _poly(c["imag"])
    return eps_real + 1j * eps_imag


def kz_soil(mv, freq_ghz=1.4, theta_inc_deg=45.0, sand=0.51, clay=0.13):
    """
    Vertical complex wavenumber in soil.

    De Zan et al. (2014) Eq. (5):
        kz' = sqrt(omega^2 * eps' * mu - kx^2)

    choosing the root with negative imaginary part (wave attenuates downward).

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture.
    freq_ghz : float
        Radar frequency in GHz.
    theta_inc_deg : float
        Incidence angle in degrees.
    sand, clay : float
        Soil texture fractions.

    Returns
    -------
    kz : complex ndarray
        Vertical wavenumber in soil [rad/m].
    """
    mv = np.asarray(mv, dtype=float)
    eps = hallikainen_dielectric(mv, freq_ghz, sand, clay)
    lam = 0.3 / freq_ghz                    # wavelength in m
    k0 = 2 * np.pi / lam                    # free-space wavenumber
    kx = k0 * np.sin(np.deg2rad(theta_inc_deg))
    kz_sq = k0**2 * eps - kx**2
    kz = np.sqrt(kz_sq.astype(complex))
    # choose root with Im(kz) < 0 (attenuation downward)
    kz = np.where(kz.imag > 0, -kz, kz)
    return kz
