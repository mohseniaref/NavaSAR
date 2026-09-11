"""
Dielectric properties of wet soil.

Hallikainen et al. (1985) empirical model for complex relative dielectric
constant as a function of volumetric soil moisture, frequency, and texture.

De Zan et al. (2014) Eq. (5): vertical complex wavenumber in soil.
"""
import numpy as np

# Hallikainen 1985, Table 1 coefficients for real (a,b,c) and imaginary (d,e,f)
# parts at 1.4 GHz (L-band) and 5.0 GHz (C-band).
# eps_r = (a0 + a1*S + a2*C) + (b0 + b1*S + b2*C)*mv + (c0 + c1*S + c2*C)*mv^2
# where S=sand fraction, C=clay fraction (0-1), mv=volumetric moisture (0-1)
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
    5.0: {
        "real": np.array([
            [2.927, -0.012, -0.001],
            [5.505,  0.371,  0.062],
            [114.826, -0.389, -0.547],
        ]),
        "imag": np.array([
            [0.004,  0.001,  0.002],
            [0.951,  0.005, -0.010],
            [15.354,  0.143, -0.261],
        ]),
    },
}


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
        Radar frequency in GHz. Supported: 1.4 (L), 5.0 (C).
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
    if freq_ghz not in _HALLIKAINEN_COEFFS:
        # linear interpolation between nearest supported frequencies
        freqs = sorted(_HALLIKAINEN_COEFFS.keys())
        f0, f1 = freqs[0], freqs[-1]
        t = (freq_ghz - f0) / (f1 - f0)
        eps0 = hallikainen_dielectric(mv, f0, sand, clay)
        eps1 = hallikainen_dielectric(mv, f1, sand, clay)
        return eps0 + t * (eps1 - eps0)

    c = _HALLIKAINEN_COEFFS[freq_ghz]
    S, C = sand, clay
    texture = np.array([1.0, S, C])

    def _poly(coeffs):
        a = coeffs @ texture          # shape (3,)
        return a[0] + a[1] * mv + a[2] * mv**2

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
