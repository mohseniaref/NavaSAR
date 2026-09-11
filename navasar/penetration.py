"""
Penetration depth and skin depth of microwaves in wet soil.

Implements:
  - De Zan (2014) Eq. (5): penetration depth from Im(kz)
  - Behari (2005) Ch. 1, Eq. (1.23): penetration depth from complex eps
  - Frequency and texture sweep helpers

References
----------
De Zan, F. et al. (2014). IEEE TGRS 52(1), 418–425.
Behari, J. (2005). Microwave Dielectric Behavior of Wet Soils. Springer.
"""
import numpy as np
from .dielectric import kz_soil, hallikainen_dielectric


# ---------------------------------------------------------------------------
# De Zan / Hallikainen penetration depth
# ---------------------------------------------------------------------------

def penetration_depth_dezan(mv, freq_ghz=1.4, theta_inc_deg=45.0,
                              sand=0.51, clay=0.13):
    """
    One-way penetration depth (1/e amplitude) from De Zan Eq. (5).

    δ = -1 / (2 Im(kz))

    where Im(kz) < 0 by convention (wave attenuates downward).

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture [m³/m³].
    freq_ghz : float
        Radar frequency in GHz.
    theta_inc_deg : float
        Incidence angle in degrees.
    sand, clay : float
        Soil texture fractions [0, 1].

    Returns
    -------
    delta : ndarray
        Penetration depth in metres.
    """
    kz = kz_soil(mv, freq_ghz, theta_inc_deg, sand, clay)
    return -1.0 / (2.0 * kz.imag)   # Im(kz) < 0


def penetration_depth_behari(eps, freq_ghz):
    """
    Penetration depth from complex permittivity.

    Behari (2005) Ch. 1, Eq. (1.23):

        Pd = λ / (4π √|ε| |sin(δ_loss)|)

    where δ_loss = angle(ε)/2 is half the loss angle.

    Parameters
    ----------
    eps : complex array_like
        Complex relative dielectric constant.
    freq_ghz : float
        Frequency in GHz.

    Returns
    -------
    Pd : ndarray
        Penetration depth in metres.
    """
    eps = np.asarray(eps, dtype=complex)
    lam = 0.3 / freq_ghz          # free-space wavelength [m]
    eps_abs = np.abs(eps)
    delta_loss = np.angle(eps) / 2.0
    Pd = lam / (4.0 * np.pi * np.sqrt(eps_abs) * np.abs(np.sin(delta_loss)))
    return np.abs(Pd)


def skin_depth_vs_moisture(mv, freq_ghz=1.4, theta_inc_deg=45.0,
                             sand=0.51, clay=0.13, method='dezan'):
    """
    Penetration depth as a function of soil moisture.

    Convenience wrapper that calls either the De Zan or Behari formula.

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture [m³/m³].
    freq_ghz : float
    theta_inc_deg : float
    sand, clay : float
    method : {'dezan', 'behari'}
        Which formula to use.

    Returns
    -------
    delta : ndarray
        Penetration depth in metres.
    """
    mv = np.asarray(mv, dtype=float)
    if method == 'dezan':
        return penetration_depth_dezan(mv, freq_ghz, theta_inc_deg, sand, clay)
    elif method == 'behari':
        eps = hallikainen_dielectric(mv, freq_ghz, sand, clay)
        return penetration_depth_behari(eps, freq_ghz)
    else:
        raise ValueError(f"method must be 'dezan' or 'behari', got '{method}'")


def penetration_depth_vs_frequency(mv, freq_array, theta_inc_deg=45.0,
                                    sand=0.51, clay=0.13):
    """
    Penetration depth as a function of frequency for fixed moisture.

    Parameters
    ----------
    mv : float
        Fixed volumetric soil moisture [m³/m³].
    freq_array : array_like
        Frequencies in GHz (e.g. np.linspace(0.5, 18, 100)).
    theta_inc_deg : float
    sand, clay : float

    Returns
    -------
    freq_array : ndarray
    delta : ndarray
        Penetration depth in metres.
    """
    freq_array = np.asarray(freq_array, dtype=float)
    delta = np.array([
        float(penetration_depth_dezan(mv, f, theta_inc_deg, sand, clay))
        for f in freq_array
    ])
    return freq_array, delta


def skin_depth_texture_map(mv, freq_ghz=1.4, theta_inc_deg=45.0,
                             n_grid=20):
    """
    2-D map of penetration depth over the sand–clay texture space.

    Parameters
    ----------
    mv : float
        Fixed volumetric soil moisture [m³/m³].
    freq_ghz : float
    theta_inc_deg : float
    n_grid : int
        Number of grid points per axis.

    Returns
    -------
    sand_grid : ndarray, shape (n_grid,)
    clay_grid : ndarray, shape (n_grid,)
    delta_map : ndarray, shape (n_grid, n_grid)
        Penetration depth [m]; NaN where sand + clay > 1.
    """
    sand_grid = np.linspace(0.05, 0.90, n_grid)
    clay_grid = np.linspace(0.05, 0.90, n_grid)
    delta_map = np.full((n_grid, n_grid), np.nan)

    for i, S in enumerate(sand_grid):
        for j, C in enumerate(clay_grid):
            if S + C > 1.0:
                continue
            delta_map[i, j] = float(
                penetration_depth_dezan(mv, freq_ghz, theta_inc_deg, S, C)
            )
    return sand_grid, clay_grid, delta_map


def two_way_penetration(mv, freq_ghz=1.4, theta_inc_deg=45.0,
                          sand=0.51, clay=0.13):
    """
    Two-way (round-trip) penetration depth for radar sensitivity estimate.

    δ_2way = δ_1way / 2   (6 dB depth)

    Parameters
    ----------
    mv : array_like
    freq_ghz, theta_inc_deg, sand, clay : float

    Returns
    -------
    delta_2way : ndarray
        Two-way penetration depth in metres.
    """
    return penetration_depth_dezan(mv, freq_ghz, theta_inc_deg, sand, clay) / 2.0
