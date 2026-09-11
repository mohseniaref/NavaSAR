"""
Dielectric models from:

  Behari, J. (2005). Microwave Dielectric Behavior of Wet Soils.
  Springer. ISBN 1-4020-3271-4.

Implements:
  - Debye equation for free water (Ch. 2, Eq. 2.12-2.14)
  - Cole-Cole equation for bound/free water (Ch. 2, Eq. 2.15)
  - Klein & Swift static dielectric constant (Ch. 2, Eq. 2.17)
  - Stogryn relaxation time (Ch. 2, Eq. 2.18)
  - Emissivity from dielectric constant (Ch. 1, Eq. 1.27)
  - Brightness temperature (Ch. 1, Eq. 1.31-1.33)
  - Wang & Schmugge mixing model (Ch. 6, Eq. 6.1-6.5)
  - Dobson semi-empirical mixing model (Ch. 6, Eq. 6.6-6.10)
  - Rayleigh mixing formula (Ch. 6, Eq. 6.24)
  - Bottcher / Polder-van Santen formula (Ch. 6, Eq. 6.33-6.34)
  - Refractive mixing formula (Ch. 6, Eq. 6.12-6.13)
  - General 4-component mixing model (Ch. 6, Eq. 6.39)
  - Penetration depth (Ch. 1, Eq. 1.23)
"""
import numpy as np


# ---------------------------------------------------------------------------
# Water dielectric models (Chapter 2)
# ---------------------------------------------------------------------------

def debye_water(freq_ghz, T_celsius=20.0):
    """
    Debye equation for free water dielectric constant.

    Behari (2005) Ch. 2, Eq. (2.12-2.14).

    Parameters
    ----------
    freq_ghz : array_like
        Frequency in GHz.
    T_celsius : float
        Temperature in °C.

    Returns
    -------
    eps : complex ndarray
        Complex relative dielectric constant of free water.
    """
    freq_ghz = np.asarray(freq_ghz, dtype=float)
    omega = 2 * np.pi * freq_ghz * 1e9

    eps_s = klein_swift_static(T_celsius)
    eps_inf = 4.2
    tau = stogryn_relaxation(T_celsius)

    eps_real = eps_inf + (eps_s - eps_inf) / (1 + omega**2 * tau**2)
    eps_imag = (eps_s - eps_inf) * omega * tau / (1 + omega**2 * tau**2)
    return eps_real + 1j * eps_imag


def cole_cole_water(freq_ghz, T_celsius=20.0, bound=False):
    """
    Cole-Cole equation for water dielectric constant.

    Behari (2005) Ch. 2, Eq. (2.15). Most accurate model for both
    free water (alpha~0.012) and bound water (alpha=0, longer tau).

    Parameters
    ----------
    freq_ghz : array_like
        Frequency in GHz.
    T_celsius : float
        Temperature in °C.
    bound : bool
        If True, use bound water parameters (Table 2.1).
        If False, use free water parameters.

    Returns
    -------
    eps : complex ndarray
        Complex relative dielectric constant.
    """
    freq_ghz = np.asarray(freq_ghz, dtype=float)
    omega = 2 * np.pi * freq_ghz * 1e9

    if bound:
        # Bound water parameters at 20°C (Behari Table 2.1)
        eps_s   = 57.9
        eps_inf = 3.15
        tau     = 9.3e-6      # much longer than free water
        alpha   = 0.0
    else:
        # Free water parameters (Behari Table 2.2 at given T)
        eps_s   = klein_swift_static(T_celsius)
        eps_inf = _eps_inf_water(T_celsius)
        tau     = stogryn_relaxation(T_celsius)
        alpha   = 0.012

    eps = eps_inf + (eps_s - eps_inf) / (1 + (1j * omega * tau)**(1 - alpha))
    # Ensure imaginary part is positive (lossy medium: eps = eps' + j*eps'')
    return eps.real + 1j * np.abs(eps.imag)


def klein_swift_static(T_celsius):
    """
    Static dielectric constant of pure water.

    Behari (2005) Ch. 2, Eq. (2.17) — Klein & Swift (1977).

    Parameters
    ----------
    T_celsius : float or array_like
        Temperature in °C.

    Returns
    -------
    eps_s : float or ndarray
    """
    T = np.asarray(T_celsius, dtype=float)
    return 88.045 - 0.4147*T + 6.295e-4*T**2 + 1.075e-5*T**3


def stogryn_relaxation(T_celsius):
    """
    Relaxation time of pure water.

    Behari (2005) Ch. 2, Eq. (2.18) — Stogryn (1970).

    Parameters
    ----------
    T_celsius : float or array_like
        Temperature in °C.

    Returns
    -------
    tau : float or ndarray
        Relaxation time in seconds.
    """
    T = np.asarray(T_celsius, dtype=float)
    two_pi_tau = (1.1109e-10 - 3.824e-12*T + 6.938e-14*T**2
                  - 5.096e-16*T**3)
    return two_pi_tau / (2 * np.pi)


def _eps_inf_water(T_celsius):
    """High-frequency limit of water dielectric constant (Table 2.2)."""
    # Interpolate from Table 2.2
    T_table = np.array([0, 10, 20, 30, 40, 50, 60])
    e_table = np.array([4.46, 4.10, 4.23, 4.20, 4.16, 4.13, 4.21])
    return float(np.interp(T_celsius, T_table, e_table))


def saline_water_dielectric(freq_ghz, T_celsius=20.0, salinity_ppt=35.0):
    """
    Complex dielectric constant of saline (sea) water.

    Behari (2005) Ch. 2, Eq. (2.19-2.21) — Klein & Swift (1977) with
    salinity correction for static constant and ionic conductivity.

    Parameters
    ----------
    freq_ghz : array_like
        Frequency in GHz.
    T_celsius : float
        Temperature in °C.
    salinity_ppt : float
        Salinity in parts per thousand (g/kg). Typical seawater ≈ 35 ppt.

    Returns
    -------
    eps : complex ndarray
        Complex relative dielectric constant of saline water.

    Notes
    -----
    Salinity reduces the static dielectric constant and adds an ionic
    conductivity term σ/(ε₀ ω) to the imaginary part, dominating at low
    frequencies (L-band).
    """
    freq_ghz = np.asarray(freq_ghz, dtype=float)
    omega = 2 * np.pi * freq_ghz * 1e9
    S_ppt = salinity_ppt           # salinity in ppt (for conductivity)
    S_frac = salinity_ppt / 1000.0 # salinity as fraction (for eps/tau formulae)
    T = T_celsius

    # Behari Eq. (2.19): salinity-corrected static dielectric constant
    # Coefficients use S as a fraction [0, 1]
    eps_s_fw = klein_swift_static(T)
    S = S_frac
    a0 = 1.0 - 0.2551 * S + 5.151e-3 * S**2 - 6.889e-5 * S**3
    eps_s = eps_s_fw * a0                                  # Eq. (2.19)

    # Behari Eq. (2.20): salinity-corrected relaxation time
    # Coefficients use S as a fraction [0, 1]
    tau_fw = stogryn_relaxation(T)
    b0 = (0.1463e-2 * S * T + 1.0
          - 0.04896 * S - 0.02967 * S**2 + 5.644e-4 * S**3)
    tau = tau_fw * b0                                      # Eq. (2.20)

    eps_inf = 4.9   # high-frequency limit

    # Debye term
    eps_debye = eps_inf + (eps_s - eps_inf) / (1.0 + 1j * omega * tau)

    # Behari Eq. (2.21): ionic conductivity sigma [S/m] — uses S in ppt
    # Stogryn (1971) conductivity formula
    S = S_ppt
    sigma = S * (0.18252 - 1.4619e-3 * S + 2.093e-5 * S**2 - 1.282e-7 * S**3)
    sigma *= np.exp(0.02033 * (T - 25) + 1.997e-4 * (T - 25)**2)
    eps0_const = 8.854e-12    # permittivity of free space [F/m]
    eps_ionic = sigma / (eps0_const * omega)

    # eps'' total = Debye imaginary + ionic conduction loss
    # eps'  total = Debye real (unchanged by ionic contribution)
    eps = eps_debye.real + 1j * (np.abs(eps_debye.imag) + eps_ionic)
    return eps.astype(complex)


# ---------------------------------------------------------------------------
# Soil physical parameters (Chapter 1)
# ---------------------------------------------------------------------------

def field_capacity(sand_pct, clay_pct):
    """
    Field capacity from soil texture.

    Behari (2005) Ch. 1, Eq. (1.8) — Schmugge et al. (1976).

    Parameters
    ----------
    sand_pct, clay_pct : float
        Sand and clay percentages (0-100).

    Returns
    -------
    FC : float
        Volumetric field capacity [m³/m³].
    """
    return (25.1 - 0.21 * sand_pct + 0.22 * clay_pct) / 100.0


def wilting_point(sand_pct, clay_pct):
    """
    Wilting point from soil texture.

    Behari (2005) Ch. 1, Eq. (1.10) — Wang & Schmugge (1980).
    """
    return 0.06774 - 0.00064 * sand_pct + 0.00478 * clay_pct


def transition_moisture(sand_pct, clay_pct):
    """
    Transition moisture (bound→free water boundary).

    Behari (2005) Ch. 1, Eq. (1.11) — Wang & Schmugge (1980).
    """
    Wp = wilting_point(sand_pct, clay_pct)
    return 0.49 * Wp + 0.165


def penetration_depth(eps, freq_ghz):
    """
    Penetration depth (1/e power) in soil.

    Behari (2005) Ch. 1, Eq. (1.23).

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
    lam = 0.3 / freq_ghz   # wavelength in m
    eps_abs = np.abs(eps)
    delta = np.angle(eps) / 2   # loss angle / 2
    Pd = lam / (4 * np.pi * np.sqrt(eps_abs) * np.abs(np.sin(delta)))
    return np.abs(Pd)


def emissivity_smooth(eps_real):
    """
    Emissivity of smooth bare soil at normal incidence.

    Behari (2005) Ch. 1, Eq. (1.27).

    Parameters
    ----------
    eps_real : array_like
        Real part of dielectric constant.

    Returns
    -------
    e : ndarray
        Emissivity [0, 1].
    """
    eps_real = np.asarray(eps_real, dtype=float)
    return 1 - ((1 - np.sqrt(eps_real)) / (1 + np.sqrt(eps_real)))**2


def emissivity_rough(eps_real, h_roughness):
    """
    Emissivity of rough soil surface.

    Behari (2005) Ch. 1, Eq. (1.30) — Choudhury et al. (1982).

    Parameters
    ----------
    eps_real : array_like
        Real part of dielectric constant.
    h_roughness : float
        Roughness parameter h (dimensionless). Smooth: h~0.1, tilled: h~0.5.

    Returns
    -------
    e_surface : ndarray
    """
    e = emissivity_smooth(eps_real)
    return 1 + (e - 1) * np.exp(h_roughness)


def brightness_temperature(eps_real, T_soil_K, h_roughness=0.0,
                            T_sky_K=5.0, T_atm_K=270.0, transmissivity=0.99):
    """
    Microwave brightness temperature of bare soil.

    Behari (2005) Ch. 1, Eq. (1.31-1.33).

    Parameters
    ----------
    eps_real : array_like
        Real part of soil dielectric constant.
    T_soil_K : float
        Soil physical temperature [K].
    h_roughness : float
        Surface roughness parameter.
    T_sky_K : float
        Sky brightness temperature [K].
    T_atm_K : float
        Average atmospheric temperature [K].
    transmissivity : float
        Atmospheric transmissivity (≈0.99 at L-band).

    Returns
    -------
    T_B : ndarray
        Brightness temperature [K].
    """
    e = emissivity_rough(eps_real, h_roughness)
    r = 1 - e   # reflectivity
    T_B = transmissivity * (r * T_sky_K + e * T_soil_K) + (1 - transmissivity) * T_atm_K
    return T_B


# ---------------------------------------------------------------------------
# Dielectric mixing models (Chapter 6)
# ---------------------------------------------------------------------------

def wang_schmugge(mv, freq_ghz=1.4, sand_pct=51.5, clay_pct=13.5,
                  rho_b=1.3, T_celsius=20.0):
    """
    Wang & Schmugge (1980) empirical mixing model.

    Behari (2005) Ch. 6, Eq. (6.1-6.5).

    Two-region model: below and above transition moisture m_t.

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture [0, 0.5].
    freq_ghz : float
        Frequency in GHz (1.4 or 5.0).
    sand_pct, clay_pct : float
        Soil texture percentages.
    rho_b : float
        Bulk density [g/cm³].
    T_celsius : float
        Temperature [°C].

    Returns
    -------
    eps : complex ndarray
        Complex relative dielectric constant of wet soil.
    """
    mv = np.asarray(mv, dtype=float)

    # Soil physical parameters
    rho_r = 2.66          # rock density g/cm³
    P = 1 - rho_b / rho_r  # porosity

    Wp = wilting_point(sand_pct, clay_pct)
    mt = transition_moisture(sand_pct, clay_pct)
    Y  = -0.57 * Wp + 0.481

    # Component dielectric constants
    eps_a = 1.0 + 0j                    # air
    eps_r = 5.5 + 0.2j                  # rock/soil solids
    eps_i = 3.2 + 0.1j                  # ice (bound water proxy)
    eps_w = cole_cole_water(freq_ghz, T_celsius, bound=False)

    # Conductivity loss term (Eq. 6.5): alpha chosen empirically
    alpha_cond = 0.3

    eps = np.zeros_like(mv, dtype=complex)

    # Region 1: mv < mt
    mask1 = mv < mt
    if np.any(mask1):
        Wc = mv[mask1]
        eps_x = eps_i + (eps_w - eps_i) * (Wc / mt) * Y
        eps[mask1] = (Wc * eps_x + (P - Wc) * eps_a + (1 - P) * eps_r)
        eps[mask1].imag += alpha_cond * Wc**2

    # Region 2: mv >= mt
    mask2 = ~mask1
    if np.any(mask2):
        Wc = mv[mask2]
        eps_x = eps_i + (eps_w - eps_i) * Y
        eps[mask2] = (mt * eps_x + (Wc - mt) * eps_w
                      + (P - Wc) * eps_a + (1 - P) * eps_r)
        eps[mask2].imag += alpha_cond * Wc**2

    return eps


def dobson_semiemp(mv, freq_ghz=1.4, sand_pct=51.5, clay_pct=13.5,
                   rho_b=1.3, T_celsius=20.0, beta=1.09, alpha=0.65):
    """
    Dobson et al. (1985) semi-empirical 4-component mixing model.

    Behari (2005) Ch. 6, Eq. (6.6-6.10).

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture [0, 0.5].
    freq_ghz : float
        Frequency in GHz.
    sand_pct, clay_pct : float
        Soil texture percentages.
    rho_b : float
        Bulk density [g/cm³].
    T_celsius : float
        Temperature [°C].
    beta : float
        Empirical parameter (~1.09 for most soils, range 1.0-1.16).
    alpha : float
        Mixing exponent (typically 0.65).

    Returns
    -------
    eps : complex ndarray
        Complex relative dielectric constant.
    """
    mv = np.asarray(mv, dtype=float)

    rho_ss = 2.65   # soil solid density g/cm³
    eps_ss = 4.7 + 0j
    eps_fw = cole_cole_water(freq_ghz, T_celsius, bound=False)

    # Eq. (6.10): eps^alpha = 1 + (rho_b/rho_ss)*(eps_ss^alpha - 1) + mv^beta*(eps_fw^alpha - 1)
    eps_alpha = (1.0
                 + (rho_b / rho_ss) * (eps_ss**alpha - 1)
                 + mv**beta * (eps_fw**alpha - 1))
    return eps_alpha**(1.0 / alpha)


def rayleigh_mixing(eps_0, eps_1, f1):
    """
    Rayleigh mixing formula for spherical inclusions.

    Behari (2005) Ch. 6, Eq. (6.24).

    Parameters
    ----------
    eps_0 : complex
        Background permittivity.
    eps_1 : complex
        Inclusion permittivity.
    f1 : array_like
        Volume fraction of inclusions [0, 1].

    Returns
    -------
    eps_eff : complex ndarray
        Effective permittivity.
    """
    f1 = np.asarray(f1, dtype=float)
    # Solve: (eps_eff - eps_0)/(eps_eff + 2*eps_0) = f1*(eps_1 - eps_0)/(eps_1 + 2*eps_0)
    rhs = f1 * (eps_1 - eps_0) / (eps_1 + 2 * eps_0)
    eps_eff = eps_0 * (1 + 2 * rhs) / (1 - rhs)
    return eps_eff


def bottcher_mixing(eps_0, eps_1, f1, tol=1e-8, max_iter=100):
    """
    Bottcher / Polder-van Santen mixing formula (v=2 case).

    Behari (2005) Ch. 6, Eq. (6.33-6.34). Implicit — solved iteratively.

    Parameters
    ----------
    eps_0 : complex
        Background permittivity.
    eps_1 : complex
        Inclusion permittivity.
    f1 : array_like
        Volume fraction of inclusions.
    tol : float
        Convergence tolerance.
    max_iter : int
        Maximum iterations.

    Returns
    -------
    eps_eff : complex ndarray
    """
    f1 = np.asarray(f1, dtype=float)
    eps_eff = np.full_like(f1, eps_0, dtype=complex)

    for _ in range(max_iter):
        eps_new = eps_0 + f1 * (eps_1 - eps_0) * 3 * eps_eff / (eps_1 + 2 * eps_eff)
        if np.max(np.abs(eps_new - eps_eff)) < tol:
            break
        eps_eff = eps_new

    return eps_eff


def refractive_mixing(eps_components, f_components):
    """
    Refractive (square-root) mixing formula for multi-component mixtures.

    Behari (2005) Ch. 6, Eq. (6.12-6.13).

    sqrt(eps_m) = sum_i f_i * sqrt(eps_i)

    Parameters
    ----------
    eps_components : list of complex
        Dielectric constants of each component.
    f_components : list of array_like
        Volume fractions of each component (must sum to 1).

    Returns
    -------
    eps_eff : complex ndarray
    """
    result = sum(np.asarray(f, dtype=float) * np.sqrt(complex(e))
                 for e, f in zip(eps_components, f_components))
    return result**2


def four_component_mixing(mv, freq_ghz=1.4, sand_pct=51.5, clay_pct=13.5,
                          rho_b=1.3, T_celsius=20.0,
                          X_bw=0.07, u=2.0, v=0.0):
    """
    General 4-component mixing model for wet soil.

    Behari (2005) Ch. 6, Eq. (6.39) — Sahu (1998), Behari et al. (2001).

    Components: dry soil (background) + free water + bound water + air.

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture [0, 0.5].
    freq_ghz : float
        Frequency in GHz.
    sand_pct, clay_pct : float
        Soil texture percentages.
    rho_b : float
        Bulk density [g/cm³].
    T_celsius : float
        Temperature [°C].
    X_bw : float
        Fraction of total water that is bound (0.05-0.10).
    u : float
        Shape factor (2 for spherical inclusions).
    v : float
        Interaction factor (0 = no neighbor influence).

    Returns
    -------
    eps_eff : complex ndarray
        Effective dielectric constant of wet soil.
    """
    mv = np.asarray(mv, dtype=float)

    rho_r = 2.66
    porosity = 1 - rho_b / rho_r

    eps_dry  = 4.7 + 0.05j                                    # Eq. (6.40)
    eps_fw   = cole_cole_water(freq_ghz, T_celsius, bound=False)
    eps_bw   = cole_cole_water(freq_ghz, T_celsius, bound=True)
    eps_a    = 1.0 + 0j

    f_bw = mv * X_bw
    f_fw = mv * (1 - X_bw)
    f_a  = np.clip(porosity - mv, 0, None)

    def _term(eps_i, f_i, eps_eff):
        num = eps_i - eps_dry
        den = (eps_eff + u * eps_dry) + v * (eps_eff - eps_dry)
        return f_i * num / den

    # Iterative solution of Eq. (6.39)
    eps_eff = eps_dry * np.ones_like(mv, dtype=complex)
    for _ in range(200):
        lhs_denom = (eps_eff + u * eps_dry) + v * (eps_eff - eps_dry)
        lhs = (eps_eff - eps_dry) / lhs_denom
        rhs = (_term(eps_fw, f_fw, eps_eff)
               + _term(eps_bw, f_bw, eps_eff)
               + _term(eps_a,  f_a,  eps_eff))
        # Update: solve for eps_eff from lhs = rhs
        # lhs ≈ rhs → eps_eff_new from rhs
        eps_new = eps_dry + rhs * ((eps_eff + u * eps_dry) + v * (eps_eff - eps_dry))
        if np.max(np.abs(eps_new - eps_eff)) < 1e-8:
            break
        eps_eff = 0.5 * eps_eff + 0.5 * eps_new   # damped update

    return eps_eff
