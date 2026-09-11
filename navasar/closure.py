"""
Closure phase time-series analysis.

Zheng & Fattahi (2026): difference between bandwidth-1 (nearest-neighbour)
and full-network InSAR time-series.
"""
import numpy as np
from .slc import simulate_slc, multilook_interferogram


def _bw1_timeseries(ifg_matrix):
    """Cumulative phase from nearest-neighbour pairs only."""
    T = ifg_matrix.shape[0]
    ts = np.zeros(T)
    for t in range(1, T):
        ts[t] = ts[t - 1] + np.angle(ifg_matrix[t - 1, t])
    return ts


def _fullnet_timeseries(ifg_matrix):
    """
    Full-network time-series via least-squares inversion.
    Solves: A * phi = dphi  (phase differences -> cumulative phases)
    """
    T = ifg_matrix.shape[0]
    n_ifg = T * (T - 1) // 2
    A = np.zeros((n_ifg, T - 1))
    dphi = np.zeros(n_ifg)
    row = 0
    for i in range(T):
        for j in range(i + 1, T):
            dphi[row] = np.angle(ifg_matrix[i, j])
            for k in range(i, j):
                A[row, k] = 1
            row += 1
    phi, _, _, _ = np.linalg.lstsq(A, dphi, rcond=None)
    return np.concatenate([[0], phi])


def closure_phase_timeseries(mv_series, freq_ghz=1.4, sand=0.51, clay=0.13,
                              n_pixels=400, n_layers=200, max_depth=0.30,
                              rng=None):
    """
    Compute closure phase time-series from a soil moisture time-series.

    Zheng & Fattahi (2026): closure_phase_ts = bw1_ts - fullnet_ts.

    Parameters
    ----------
    mv_series : array_like, shape (T,) or (T, N)
        Volumetric soil moisture per acquisition [and per layer].
    freq_ghz : float
    sand, clay : float
    n_pixels : int
    n_layers : int
    max_depth : float
    rng : np.random.Generator or None

    Returns
    -------
    closure_ts : ndarray, shape (T,)
        Closure phase time-series in degrees.
    bw1_ts : ndarray
    fullnet_ts : ndarray
    """
    mv_or_slc = np.asarray(mv_series)
    if np.iscomplexobj(mv_or_slc) and mv_or_slc.ndim == 2:
        slc_stack = mv_or_slc
    else:
        slc_stack = simulate_slc(mv_series, freq_ghz=freq_ghz, sand=sand,
                                  clay=clay, n_pixels=n_pixels, n_layers=n_layers,
                                  max_depth=max_depth, rng=rng)
    T = slc_stack.shape[0]

    # Build interferogram matrix
    ifg_matrix = np.zeros((T, T), dtype=complex)
    for i in range(T):
        for j in range(i + 1, T):
            ifg_matrix[i, j] = multilook_interferogram(slc_stack[i], slc_stack[j])
            ifg_matrix[j, i] = np.conj(ifg_matrix[i, j])

    bw1 = _bw1_timeseries(ifg_matrix)
    full = _fullnet_timeseries(ifg_matrix)

    closure_ts = np.rad2deg(bw1 - full)
    return closure_ts, np.rad2deg(bw1), np.rad2deg(full)
