"""
Single-look complex (SLC) SAR simulation.

Zheng & Fattahi (2026) Eq. (3):
    s1 = sum_n theta_n * exp(-2j * k1n * n * D/N)

where theta_n ~ CN(0, sigma_n).
"""
import numpy as np
from .dielectric import kz_soil


def _make_sigma(n_layers, max_depth, profile, scale_depth=0.10):
    """
    Scatterer strength (sigma) profile over depth.

    Parameters
    ----------
    n_layers : int
    max_depth : float
        Total depth [m].
    profile : {'uniform', 'exponential'}
        'uniform'     — equal weight at all depths (Zheng default).
        'exponential' — sigma_n ∝ exp(-depth / scale_depth), mimicking
                        rapid attenuation; highlights near-surface layers.
    scale_depth : float
        e-folding depth [m] for the exponential profile.

    Returns
    -------
    sigma : ndarray, shape (n_layers,)
    """
    depths = np.arange(n_layers) * max_depth / n_layers
    if profile == 'uniform':
        sigma = np.ones(n_layers)
    elif profile == 'exponential':
        sigma = np.exp(-depths / scale_depth)
    else:
        raise ValueError(f"sigma_profile must be 'uniform' or 'exponential', "
                         f"got '{profile}'")
    return sigma / sigma.sum()   # normalise so total power is constant


def simulate_slc(mv_series, freq_ghz=1.4, sand=0.51, clay=0.13,
                 n_pixels=400, n_layers=200, max_depth=0.30,
                 sigma=None, sigma_profile='uniform', scale_depth=0.10,
                 rng=None):
    """
    Simulate a stack of SLC SAR acquisitions.

    Zheng & Fattahi (2026) Eq. (3), Table 1.

    Parameters
    ----------
    mv_series : array_like, shape (T,) or (T, N)
        Volumetric soil moisture per acquisition [and per layer].
        If 1-D, uniform moisture with depth is assumed.
    freq_ghz : float
    sand, clay : float
    n_pixels : int
        Number of pixels in the multi-look window (M).
    n_layers : int
        Number of depth layers (N).
    max_depth : float
        Maximum penetration depth in metres (D).
    sigma : array_like or None
        Explicit radar cross-section profile per layer.
        If None, ``sigma_profile`` is used to generate it.
    sigma_profile : {'uniform', 'exponential'}
        Scatterer strength profile when ``sigma`` is None.
        'uniform'     — all layers equally weighted (Zheng & Fattahi 2026
                        default, Table 1).
        'exponential' — strength decays as exp(-depth / scale_depth),
                        emphasising near-surface scattering (Fig. 7 variant).
    scale_depth : float
        e-folding depth [m] for the exponential profile. Default 0.10 m.
    rng : np.random.Generator or None

    Returns
    -------
    slc_stack : complex ndarray, shape (T, n_pixels)
    """
    if rng is None:
        rng = np.random.default_rng(42)

    mv_series = np.asarray(mv_series, dtype=float)
    T = mv_series.shape[0]

    if mv_series.ndim == 1:
        mv_series = np.tile(mv_series[:, None], (1, n_layers))
    else:
        n_layers = mv_series.shape[1]

    depths = np.arange(n_layers) * max_depth / n_layers   # shape (N,)

    if sigma is None:
        sigma = _make_sigma(n_layers, max_depth, sigma_profile, scale_depth)

    # theta_n ~ CN(0, sigma_n): shape (n_pixels, n_layers)
    theta = (rng.standard_normal((n_pixels, n_layers)) +
             1j * rng.standard_normal((n_pixels, n_layers))) / np.sqrt(2)
    theta *= np.sqrt(sigma)[None, :]

    slc_stack = np.zeros((T, n_pixels), dtype=complex)
    for t in range(T):
        kz = kz_soil(mv_series[t], freq_ghz, sand=sand, clay=clay)  # shape (N,)
        phase = np.exp(-2j * kz[None, :] * depths[None, :])          # (1, N)
        slc_stack[t] = (theta * phase).sum(axis=1)                    # (n_pixels,)

    return slc_stack


def multilook_interferogram(slc1, slc2):
    """
    Multi-looked interferogram from two SLC arrays.

    Zheng & Fattahi (2026) Eq. (4).

    Parameters
    ----------
    slc1, slc2 : complex ndarray, shape (n_pixels,)

    Returns
    -------
    z12 : complex scalar
    """
    return np.mean(slc1 * np.conj(slc2))
