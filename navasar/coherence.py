"""
Interferometric coherence and phase triplets.

De Zan et al. (2014) Eq. (12-14): complex coherence for uniform vertical
scattering profile.

Phase triplet (closure phase) definition: Eq. (15).
"""
import numpy as np
from .dielectric import kz_soil


def interferometric_coherence(mv1, mv2, freq_ghz=1.4, theta_inc_deg=45.0,
                               sand=0.51, clay=0.13):
    """
    Complex interferometric coherence between two acquisitions.

    De Zan et al. (2014) Eq. (14):
        gamma = (2j * Im(kz1) * Im(kz2))^0.5 / (kz1* - 2*kz2*)

    Assumes uniform vertical scattering profile (Eq. 12).

    Parameters
    ----------
    mv1, mv2 : array_like
        Volumetric soil moisture of master and slave.
    freq_ghz : float
    theta_inc_deg : float
    sand, clay : float

    Returns
    -------
    gamma : complex ndarray
        Complex coherence (magnitude <= 1).
    """
    mv1 = np.asarray(mv1, dtype=float)
    mv2 = np.asarray(mv2, dtype=float)
    kz1 = kz_soil(mv1, freq_ghz, theta_inc_deg, sand, clay)
    kz2 = kz_soil(mv2, freq_ghz, theta_inc_deg, sand, clay)

    # Eq. (12): I(eps1, eps2) = 1 / (2j*kz1 - 2j*kz2*)
    I12 = 1.0 / (2j * kz1 - 2j * np.conj(kz2))
    I11 = 1.0 / (2j * kz1 - 2j * np.conj(kz1))
    I22 = 1.0 / (2j * kz2 - 2j * np.conj(kz2))

    gamma = I12 / np.sqrt(I11 * I22)
    return gamma


def phase_triplet(mv1, mv2, mv3, freq_ghz=1.4, theta_inc_deg=45.0,
                  sand=0.51, clay=0.13):
    """
    Interferometric phase triplet (closure phase).

    De Zan et al. (2014) Eq. (15):
        phi_123 = phi_12 + phi_23 - phi_13

    Parameters
    ----------
    mv1, mv2, mv3 : array_like
        Soil moisture at three acquisition times.

    Returns
    -------
    triplet : ndarray
        Phase mismatch in degrees.
    """
    kw = dict(freq_ghz=freq_ghz, theta_inc_deg=theta_inc_deg,
              sand=sand, clay=clay)
    g12 = interferometric_coherence(mv1, mv2, **kw)
    g23 = interferometric_coherence(mv2, mv3, **kw)
    g13 = interferometric_coherence(mv1, mv3, **kw)

    # Triplet = angle of the triple product (De Zan 2014, Eq. 15)
    triplet_rad = np.angle(g12 * g23 * np.conj(g13))
    return np.rad2deg(triplet_rad)


def fresnel_phase(mv, freq_ghz=1.4, theta_inc_deg=45.0, sand=0.51, clay=0.13):
    """
    One-way phase of Fresnel transmission coefficient (TE/HH case).

    De Zan et al. (2014) Eq. (3), Fig. 5.

    Returns
    -------
    phase_te, phase_tm : ndarray in degrees
    """
    from .dielectric import hallikainen_dielectric
    mv = np.asarray(mv, dtype=float)
    eps = hallikainen_dielectric(mv, freq_ghz, sand, clay)
    lam = 0.3 / freq_ghz
    k0 = 2 * np.pi / lam
    theta = np.deg2rad(theta_inc_deg)
    kz_air = k0 * np.cos(theta)
    kx = k0 * np.sin(theta)
    kz_soil_val = np.sqrt((k0**2 * eps - kx**2).astype(complex))
    kz_soil_val = np.where(kz_soil_val.imag > 0, -kz_soil_val, kz_soil_val)

    tau_te = 2 * kz_air / (kz_air + kz_soil_val)
    tau_tm = 2 * eps * kz_air / (eps * kz_air + kz_soil_val)
    # Return phase relative to dry soil (mv=0) so the plot shows variation only
    eps0   = hallikainen_dielectric(np.zeros_like(mv), freq_ghz, sand, clay)
    kz0    = np.sqrt((k0**2 * eps0 - kx**2).astype(complex))
    kz0    = np.where(kz0.imag > 0, -kz0, kz0)
    tau_te0 = 2 * kz_air / (kz_air + kz0)
    tau_tm0 = 2 * eps0 * kz_air / (eps0 * kz_air + kz0)
    return (np.rad2deg(np.angle(tau_te)) - np.rad2deg(np.angle(tau_te0)),
            np.rad2deg(np.angle(tau_tm)) - np.rad2deg(np.angle(tau_tm0)))
