"""
Soil moisture retrieval from InSAR observables.

Implements:
  - De Zan et al. (2014) Sec. III-D: numerical inversion of |γ| → mv
  - Zheng & Fattahi (2026) Sec. 4: closure phase step → moisture anomaly
  - InSAR Soil Moisture Index (SMI) time-series from closure phase gradient

References
----------
De Zan, F. et al. (2014). IEEE TGRS 52(1), 418–425.
Zheng, Y. & Fattahi, H. (2026). RSE 333, 115104.
"""
import numpy as np
from scipy.optimize import brentq, least_squares
from .coherence import interferometric_coherence
from .closure import closure_phase_timeseries


# ---------------------------------------------------------------------------
# De Zan (2014) coherence-based inversion
# ---------------------------------------------------------------------------

def calibration_coherence_curve(mv_ref, freq_ghz=1.4, theta_inc_deg=45.0,
                                  sand=0.51, clay=0.13,
                                  n_points=100):
    """
    Forward model: |γ| as a function of slave moisture, master fixed at mv_ref.

    De Zan et al. (2014) Sec. III-D. Used to build a lookup table for
    numerical inversion.

    Parameters
    ----------
    mv_ref : float
        Reference (master) volumetric soil moisture [m³/m³].
    freq_ghz : float
        Radar frequency in GHz.
    theta_inc_deg : float
        Incidence angle in degrees.
    sand, clay : float
        Soil texture fractions.
    n_points : int
        Number of points in the lookup table.

    Returns
    -------
    mv_slave : ndarray, shape (n_points,)
        Slave moisture values [m³/m³].
    gamma_mag : ndarray, shape (n_points,)
        Corresponding |γ| values.
    gamma_phase_deg : ndarray, shape (n_points,)
        Corresponding phase of γ in degrees.
    """
    mv_slave = np.linspace(0.0, 0.50, n_points)
    gamma = interferometric_coherence(
        mv_ref * np.ones(n_points), mv_slave,
        freq_ghz=freq_ghz, theta_inc_deg=theta_inc_deg,
        sand=sand, clay=clay
    )
    return mv_slave, np.abs(gamma), np.rad2deg(np.angle(gamma))


def invert_moisture_from_coherence(gamma_obs, mv_ref,
                                    freq_ghz=1.4, theta_inc_deg=45.0,
                                    sand=0.51, clay=0.13,
                                    observable='magnitude'):
    """
    Retrieve slave soil moisture from observed interferometric coherence.

    De Zan et al. (2014) Sec. III-D: the forward model |γ|(mv_slave) is
    monotonically decreasing in most practical ranges; root-finding via
    Brent's method.

    Parameters
    ----------
    gamma_obs : float or array_like
        Observed complex coherence, |γ| (if observable='magnitude')
        or full complex γ (if observable='phase').
    mv_ref : float
        Reference (master) moisture [m³/m³].
    freq_ghz : float
    theta_inc_deg : float
    sand, clay : float
    observable : {'magnitude', 'phase'}
        Which component to invert. 'magnitude' is better constrained.

    Returns
    -------
    mv_retrieved : float or ndarray
        Retrieved slave moisture [m³/m³].
    residual : float or ndarray
        Residual: |model - observed| at the solution.

    Raises
    ------
    ValueError
        If inversion fails to converge (observable out of range).
    """
    gamma_obs = np.atleast_1d(np.asarray(gamma_obs, dtype=float))
    mv_retrieved = np.zeros_like(gamma_obs)
    residual     = np.zeros_like(gamma_obs)

    for i, g_obs in enumerate(gamma_obs.flat):
        def residual_fn(mv2):
            gamma_model = interferometric_coherence(
                mv_ref, float(mv2),
                freq_ghz=freq_ghz, theta_inc_deg=theta_inc_deg,
                sand=sand, clay=clay
            )
            if observable == 'magnitude':
                return float(np.abs(gamma_model)) - g_obs
            else:
                return float(np.rad2deg(np.angle(gamma_model))) - g_obs

        try:
            sol = brentq(residual_fn, 0.001, 0.499, xtol=1e-5, full_output=True)
            mv_retrieved.flat[i] = sol[0]
            residual.flat[i]     = abs(residual_fn(sol[0]))
        except ValueError:
            mv_retrieved.flat[i] = np.nan
            residual.flat[i]     = np.nan

    if mv_retrieved.size == 1:
        return float(mv_retrieved[0]), float(residual[0])
    return mv_retrieved, residual


# ---------------------------------------------------------------------------
# Zheng & Fattahi (2026) closure-phase retrieval
# ---------------------------------------------------------------------------

def calibrate_closure_transfer(mv_base, freq_ghz=1.4, sand=0.51, clay=0.13,
                                 n_pixels=400, n_layers=200, max_depth=0.30,
                                 anomaly_amps=None, tau_acq=4, rng=None):
    """
    Build a calibration curve: moisture anomaly α → closure phase step β.

    Zheng & Fattahi (2026) Eq. (8):  β ≈ T_coeff × α

    The transfer coefficient T_coeff [deg / (m³/m³)] depends on:
    frequency, soil texture, baseline moisture, and simulation parameters.

    Parameters
    ----------
    mv_base : float
        Baseline (dry) soil moisture [m³/m³].
    freq_ghz : float
    sand, clay : float
    n_pixels, n_layers : int
    max_depth : float
    anomaly_amps : array_like or None
        Amplitude values to scan. Defaults to np.linspace(0.01, 0.30, 20).
    tau_acq : int
        Duration of the moisture anomaly in acquisitions.
    rng : np.random.Generator or None

    Returns
    -------
    anomaly_amps : ndarray
        Scanned anomaly amplitudes [m³/m³].
    cp_steps : ndarray
        Corresponding closure phase steps [deg].
    T_coeff : float
        Linear transfer coefficient [deg / (m³/m³)].
    """
    if rng is None:
        rng = np.random.default_rng(0)
    if anomaly_amps is None:
        anomaly_amps = np.linspace(0.01, 0.30, 20)
    anomaly_amps = np.asarray(anomaly_amps, dtype=float)

    T_total = 30   # total acquisitions
    t_start = 8
    cp_steps = []

    for amp in anomaly_amps:
        mv_ts = np.full(T_total, mv_base)
        mv_ts[t_start:t_start + tau_acq] = np.clip(mv_base + amp, 0, 0.50)
        cp, _, _ = closure_phase_timeseries(
            mv_ts, freq_ghz=freq_ghz, sand=sand, clay=clay,
            n_pixels=n_pixels, n_layers=n_layers, max_depth=max_depth,
            rng=np.random.default_rng(rng.integers(0, 2**31))
        )
        cp_steps.append(np.max(cp) - np.min(cp))

    cp_steps = np.array(cp_steps)
    # Linear fit through origin: β = T_coeff * α
    T_coeff = float(np.polyfit(anomaly_amps, cp_steps, 1)[0])
    return anomaly_amps, cp_steps, T_coeff


def invert_moisture_from_triplet(triplet_step_deg, mv_base,
                                   freq_ghz=1.4, sand=0.51, clay=0.13,
                                   n_pixels=400, n_layers=200,
                                   max_depth=0.30, rng=None):
    """
    Retrieve moisture anomaly magnitude from observed closure phase step.

    Zheng & Fattahi (2026) Sec. 4, inverts: α = β / T_coeff

    Parameters
    ----------
    triplet_step_deg : float or array_like
        Observed closure phase step [deg].
    mv_base : float
        Baseline moisture [m³/m³].
    freq_ghz, sand, clay, n_pixels, n_layers, max_depth
        Simulation parameters (same as used to compute the forward model).
    rng : np.random.Generator or None

    Returns
    -------
    mv_anomaly : float or ndarray
        Retrieved moisture anomaly [m³/m³].
    T_coeff : float
        Transfer coefficient used [deg / (m³/m³)].
    """
    _, _, T_coeff = calibrate_closure_transfer(
        mv_base, freq_ghz=freq_ghz, sand=sand, clay=clay,
        n_pixels=n_pixels, n_layers=n_layers, max_depth=max_depth, rng=rng
    )
    triplet_step_deg = np.asarray(triplet_step_deg, dtype=float)
    mv_anomaly = triplet_step_deg / T_coeff
    return mv_anomaly, T_coeff


def insar_smi_timeseries(cp_ts_deg, T_coeff, smooth_window=3):
    """
    InSAR Soil Moisture Index (SMI) from closure phase time-series.

    Zheng & Fattahi (2026) Sec. 4: the gradient of the closure phase
    time-series is proportional to instantaneous moisture anomalies.

    Parameters
    ----------
    cp_ts_deg : array_like, shape (T,)
        Closure phase time-series [deg].
    T_coeff : float
        Transfer coefficient [deg / (m³/m³)], from calibrate_closure_transfer.
    smooth_window : int
        Uniform smoothing window length for gradient estimation (odd).

    Returns
    -------
    smi : ndarray, shape (T,)
        Soil Moisture Index [m³/m³] — relative to the first acquisition.

    Notes
    -----
    The SMI represents the instantaneous moisture anomaly relative to the
    system's background; absolute moisture requires external calibration.
    """
    cp = np.asarray(cp_ts_deg, dtype=float)
    # Smooth before differentiating to reduce noise
    w = smooth_window
    kernel = np.ones(w) / w
    cp_smooth = np.convolve(cp, kernel, mode='same')
    gradient = np.gradient(cp_smooth)
    return gradient / T_coeff


# ---------------------------------------------------------------------------
# Multi-observable joint inversion
# ---------------------------------------------------------------------------

def joint_inversion(gamma_mag_obs, triplet_obs_deg,
                     mv_ref, freq_ghz=1.4, theta_inc_deg=45.0,
                     sand=0.51, clay=0.13,
                     weights=(1.0, 1.0)):
    """
    Joint inversion of |γ| and closure phase for slave moisture.

    Minimises the weighted residual:

        R(mv2) = w1 * (|γ_model| - |γ_obs|)² + w2 * (Φ_model - Φ_obs)²

    Parameters
    ----------
    gamma_mag_obs : float
        Observed coherence magnitude |γ|.
    triplet_obs_deg : float
        Observed phase triplet or closure phase step [deg].
    mv_ref : float
        Master moisture [m³/m³].
    freq_ghz, theta_inc_deg, sand, clay : float
        Observation geometry.
    weights : tuple of float
        Relative weights (w_coherence, w_triplet).

    Returns
    -------
    mv_retrieved : float
        Best-fit slave moisture [m³/m³].
    residual : float
        Weighted residual at the solution.
    """
    w1, w2 = weights

    def cost(mv2):
        mv2 = float(mv2[0])
        gamma = interferometric_coherence(
            mv_ref, mv2,
            freq_ghz=freq_ghz, theta_inc_deg=theta_inc_deg,
            sand=sand, clay=clay
        )
        r1 = w1 * (np.abs(gamma) - gamma_mag_obs)**2
        r2 = w2 * (np.rad2deg(np.angle(gamma)) - triplet_obs_deg)**2
        return r1 + r2

    result = least_squares(cost, x0=[mv_ref], bounds=([0.0], [0.5]),
                           method='trf', xtol=1e-6)
    return float(result.x[0]), float(result.cost)
