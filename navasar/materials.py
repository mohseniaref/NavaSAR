"""Microwave dielectric building blocks based on Long & Ulaby, Chapter 4.

The module complements the soil-focused :mod:`navasar.dielectric` module.
Complex permittivity follows NavaSAR's ``eps' + 1j*eps''`` loss convention.
"""
import numpy as np

from .behari import debye_water, saline_water_dielectric

_EPS0 = 8.854e-12  # vacuum permittivity, F/m (Sec. 4-8.2)


def pure_water_permittivity(freq_ghz, temperature_c=20.0):
    """Single-Debye pure-water permittivity (Chapter 4, Sec. 4-1)."""
    return debye_water(freq_ghz, T_celsius=temperature_c)


def saline_water_permittivity(freq_ghz, temperature_c=20.0, salinity_psu=35.0):
    """Saline-water permittivity including ionic loss (Sec. 4-2)."""
    return saline_water_dielectric(freq_ghz, temperature_c, salinity_psu)


def pure_ice_permittivity(freq_ghz, temperature_c=-10.0):
    """Pure-ice permittivity, Eqs. (4.23)--(4.25).

    Valid for microwave frequencies and approximately -40 to 0 degrees C.
    Implements the Matzler/Hufford loss model with frequency in GHz.
    """
    f = np.asarray(freq_ghz, dtype=float)
    temperature = np.asarray(temperature_c, dtype=float)
    if np.any(f <= 0):
        raise ValueError("frequency must be positive")
    if np.any((temperature < -40) | (temperature > 0)):
        raise ValueError("pure-ice model is valid from -40 to 0 degrees C")
    tk = temperature + 273.15
    theta = 300.0 / tk - 1.0
    alpha0 = (0.00504 + 0.0062 * theta) * np.exp(-22.1 * theta)
    b1, b, b2 = 0.0207, 335.0, 1.16e-11
    exp_term = np.exp(b / tk)
    beta0 = (b1 / tk * exp_term / (exp_term - 1.0) ** 2
             + b2 * f ** 2
             + np.exp(-9.963 + 0.0372 * (tk - 273.16)))
    eps_real = 3.1884 + 9.1e-4 * temperature
    return eps_real + 1j * (alpha0 / f + beta0 * f)


def brine_salinity(temperature_c):
    """Equilibrium sea-ice brine salinity in psu, Eq. (4.46)."""
    t = np.asarray(temperature_c, dtype=float)
    if np.any((t > -2.0) | (t < -43.2)):
        raise ValueError("brine-salinity model is valid from -43.2 to -2 degrees C")
    return np.select(
        [t >= -8.2, t >= -22.9, t >= -36.8],
        [1.725 - 18.756*t - 0.3964*t**2,
         57.041 - 9.929*t - 0.16204*t**2 - 0.002396*t**3,
         242.94 + 1.5299*t + 0.0429*t**2],
        default=508.18 + 14.535*t + 0.2018*t**2,
    )


def brine_volume_fraction(ice_salinity_psu, temperature_c):
    """Liquid-brine volume fraction in sea ice, Eq. (4.51)."""
    t = np.asarray(temperature_c, dtype=float)
    salinity = np.asarray(ice_salinity_psu, dtype=float)
    if np.any((t > -0.5) | (t < -22.9)) or np.any(salinity < 0):
        raise ValueError("requires salinity >= 0 and temperature from -22.9 to -0.5 C")
    return 1e-3 * salinity * (-49.185 / t + 0.532)


def sea_ice_permittivity(freq_ghz, temperature_c=-10.0,
                         ice_salinity_psu=5.0, inclusion_shape="spheres"):
    """Sea-ice mixture of pure ice and equilibrium liquid brine.

    Brine volume follows Eq. (4.51); ``spheres`` uses the TVB formula and
    ``random`` uses symmetric Polder--van Santen mixing. Air pockets are not
    included, so this is a compact two-phase model rather than a full ice-type
    classification model.
    """
    vb = brine_volume_fraction(ice_salinity_psu, temperature_c)
    eps_ice = pure_ice_permittivity(freq_ghz, temperature_c)
    eps_brine = saline_water_permittivity(
        freq_ghz, temperature_c, brine_salinity(np.minimum(temperature_c, -2.0))
    )
    if inclusion_shape == "spheres":
        return tvb_mixing(eps_ice, eps_brine, vb, shape="spheres")
    if inclusion_shape == "random":
        return polder_van_santen([eps_ice, eps_brine], [1.0-vb, vb])
    raise ValueError("inclusion_shape must be 'spheres' or 'random'")


def depolarization_factors(axis_a, axis_b, axis_c, quadrature_order=256):
    """Depolarization factors of an ellipsoid (Eq. 4.27).

    The factors are evaluated by Gauss-Legendre quadrature after mapping the
    semi-infinite integration interval to [0, 1]. They sum to one.
    """
    axes = np.asarray([axis_a, axis_b, axis_c], dtype=float)
    if np.any(axes <= 0):
        raise ValueError("ellipsoid semi-axes must be positive")
    nodes, weights = np.polynomial.legendre.leggauss(int(quadrature_order))
    t = (nodes + 1.0) / 2.0
    w = weights / 2.0
    # s=t/(1-t), ds=dt/(1-t)^2
    s = t / (1.0 - t)
    jac = 1.0 / (1.0 - t) ** 2
    delta = np.sqrt(np.prod(s[:, None] + axes[None, :] ** 2, axis=1))
    prefactor = np.prod(axes) / 2.0
    factors = np.array([
        prefactor * np.sum(w * jac / ((s + axis ** 2) * delta))
        for axis in axes
    ])
    return factors / factors.sum()


def dilute_ellipsoid_mixing(eps_host, eps_inclusion, volume_fraction,
                            depolarization=(1/3, 1/3, 1/3)):
    """Randomly oriented dilute ellipsoidal inclusions (Sec. 4-4.1).

    For spherical inclusions this reduces to the Maxwell-Garnett dilute form.
    """
    f = np.asarray(volume_fraction, dtype=float)
    if np.any((f < 0) | (f > 1)):
        raise ValueError("volume_fraction must lie in [0, 1]")
    factors = np.asarray(depolarization, dtype=float)
    if factors.shape != (3,) or np.any(factors < 0) or not np.isclose(factors.sum(), 1):
        raise ValueError("depolarization must contain three non-negative factors summing to one")
    host, inclusion = complex(eps_host), complex(eps_inclusion)
    contrast = inclusion - host
    polarizability = np.mean(host * contrast / (host + factors * contrast))
    return host + f * polarizability


def polder_van_santen(eps_components, volume_fractions, tol=1e-10,
                      max_iter=500):
    """Symmetric Polder-van Santen/Bruggeman mixture for spheres (Sec. 4-4.2).

    Solves ``sum(f_i (eps_i-eps)/(eps_i+2 eps)) = 0`` by damped Newton
    iteration. Fractions must sum to one.
    """
    eps_i = np.asarray(eps_components, dtype=complex)
    fractions = np.asarray(volume_fractions, dtype=float)
    if eps_i.ndim != 1 or fractions.shape != eps_i.shape:
        raise ValueError("components and fractions must be equal-length 1-D arrays")
    if np.any(fractions < 0) or not np.isclose(fractions.sum(), 1.0):
        raise ValueError("volume fractions must be non-negative and sum to one")
    effective = np.sum(fractions * eps_i)
    initial = effective
    for _ in range(int(max_iter)):
        denom = eps_i + 2.0 * effective
        residual = np.sum(fractions * (eps_i - effective) / denom)
        derivative = np.sum(fractions * (-3.0 * eps_i) / denom ** 2)
        if abs(derivative) < 1e-15:
            break
        candidate = effective - residual / derivative
        updated = 0.5 * effective + 0.5 * candidate
        if abs(updated - effective) <= tol * max(1.0, abs(updated)):
            return updated
        effective = updated

    # High dielectric contrast can make the complex Newton path oscillatory.
    # Fall back to a bounded real/imaginary least-squares solve.
    from scipy.optimize import least_squares

    def residual_2d(value):
        trial = value[0] + 1j * value[1]
        r = np.sum(fractions * (eps_i - trial) / (eps_i + 2.0*trial))
        return [r.real, r.imag]

    scale = max(1.0, float(np.max(np.abs(eps_i))))
    solution = least_squares(
        residual_2d, [initial.real, initial.imag],
        bounds=([-scale, 0.0], [10.0*scale, 10.0*scale]),
        xtol=tol, ftol=tol, gtol=tol, max_nfev=5*int(max_iter),
    )
    result = solution.x[0] + 1j*solution.x[1]
    if solution.success and np.linalg.norm(residual_2d(solution.x)) < 100*tol:
        return result
    raise RuntimeError("Polder-van Santen solver did not converge")


def dry_snow_permittivity(density_g_cm3, ice_loss=0.0):
    """Common dry-snow density relation summarized in Sec. 4-6.

    ``eps' = 1 + 1.7 rho + 0.7 rho^2`` for density ``rho`` in g/cm^3.
    The optional loss term permits coupling to a frequency-dependent ice model.
    """
    rho = np.asarray(density_g_cm3, dtype=float)
    if np.any((rho < 0) | (rho > 0.917)):
        raise ValueError("snow density must lie between 0 and pure-ice density")
    return 1.0 + 1.7 * rho + 0.7 * rho ** 2 + 1j * np.asarray(ice_loss, dtype=float)


def dry_snow_tvb(freq_ghz, density_g_cm3, temperature_c=-10.0):
    """Dry-snow TVB model with spherical ice inclusions, Eqs. (4.53)-(4.57)."""
    rho = np.asarray(density_g_cm3, dtype=float)
    if np.any((rho < 0) | (rho > 0.917)):
        raise ValueError("snow density must lie between 0 and 0.917 g/cm3")
    ice_fraction = rho / 0.917
    return tvb_mixing(1.0, pure_ice_permittivity(freq_ghz, temperature_c),
                      ice_fraction, shape="spheres")


def wet_snow_permittivity(freq_ghz, dry_density_g_cm3,
                          liquid_water_fraction, temperature_c=0.0):
    """Three-phase wet-snow mixture (air, ice, liquid water).

    Liquid water fraction is volumetric. The remaining solid volume is fixed
    by dry-snow density and mixed with a symmetric Polder--van Santen model.
    """
    rho = float(dry_density_g_cm3)
    water = float(liquid_water_fraction)
    ice = rho / 0.917
    air = 1.0 - ice - water
    if min(air, ice, water) < 0:
        raise ValueError("air, ice, and water fractions must be non-negative")
    eps_ice = pure_ice_permittivity(freq_ghz, min(float(temperature_c), 0.0))
    eps_water = pure_water_permittivity(freq_ghz, temperature_c)
    return polder_van_santen([1.0, eps_ice, eps_water], [air, ice, water])


def powdered_rock_permittivity(solid_permittivity, bulk_density_g_cm3,
                               solid_density_g_cm3=2.65):
    """Air/rock powder mixture parameterized by bulk density (Sec. 4-7)."""
    fraction = np.asarray(bulk_density_g_cm3, dtype=float) / float(solid_density_g_cm3)
    if np.any((fraction < 0) | (fraction > 1)):
        raise ValueError("bulk density must lie between zero and solid density")
    return tvb_mixing(1.0, solid_permittivity, fraction, shape="spheres")


def bound_water_permittivity(freq_ghz):
    """Vegetation bound-water Cole-Cole model, Eq. (4.76)."""
    f = np.asarray(freq_ghz, dtype=float)
    if np.any(f <= 0):
        raise ValueError("frequency must be positive")
    raw = 2.9 + 55.0 / (1.0 + (1j * f / 0.18) ** 0.5)
    return raw.real + 1j * np.abs(raw.imag)


def vegetation_permittivity(freq_ghz, gravimetric_moisture,
                            temperature_c=22.0, salinity_psu=0.0):
    """Ulaby--El-Rayes vegetation-material model, Eqs. (4.72)-(4.78).

    Gravimetric moisture ``mg`` is a fraction in [0, 1]. The result describes
    plant material, not a canopy-scale effective medium containing air.
    """
    mg = np.asarray(gravimetric_moisture, dtype=float)
    if np.any((mg < 0) | (mg > 1)):
        raise ValueError("gravimetric moisture must lie in [0, 1]")
    residual = 1.7 - 0.74*mg + 6.16*mg**2
    free_fraction = np.maximum(0.0, mg * (0.55*mg - 0.076))
    bound_fraction = 4.64*mg**2 / (1.0 + 7.36*mg**2)
    free = saline_water_permittivity(freq_ghz, temperature_c, salinity_psu)
    bound = bound_water_permittivity(freq_ghz)
    return residual + free_fraction*free + bound_fraction*bound


def canopy_effective_permittivity(plant_permittivity, vegetation_fraction,
                                  shape="needles"):
    """Mix vegetation material into air to obtain canopy-scale permittivity."""
    return tvb_mixing(1.0, plant_permittivity, vegetation_fraction, shape=shape)


def refractive_mixture(eps_components, volume_fractions):
    """Refractive-index mixture ``sqrt(eps)=sum(f_i sqrt(eps_i))``."""
    eps_i = np.asarray(eps_components, dtype=complex)
    fractions = np.asarray(volume_fractions, dtype=float)
    if fractions.shape != eps_i.shape or np.any(fractions < 0):
        raise ValueError("components and non-negative fractions must have matching shapes")
    if not np.isclose(fractions.sum(), 1.0):
        raise ValueError("volume fractions must sum to one")
    return np.sum(fractions * np.sqrt(eps_i)) ** 2


def tvb_mixing(eps_host, eps_inclusion, volume_fraction, shape="spheres"):
    """Tinga--Voss--Blossey two-phase formulas, Eqs. (4.42)-(4.44)."""
    host, inclusion = np.asarray(eps_host, dtype=complex), np.asarray(eps_inclusion, dtype=complex)
    f = np.asarray(volume_fraction, dtype=float)
    if np.any((f < 0) | (f > 1)):
        raise ValueError("volume_fraction must lie in [0, 1]")
    if shape == "spheres":
        return host + 3*f*host*(inclusion-host) / ((2*host+inclusion)-f*(inclusion-host))
    if shape == "discs":
        factor = (2*inclusion*(1-f)+host*(1+2*f)) / (f*host+(1-f)*inclusion)
        return host + f*(inclusion-host)*factor/3
    if shape == "needles":
        factor = (host*(5+f)+(1-f)*inclusion) / (host*(1+f)+inclusion*(1-f))
        return host + f*(inclusion-host)*factor/3
    raise ValueError("shape must be 'spheres', 'discs', or 'needles'")


def generalized_dielectric_mixing(eps_host, eps_inclusion, volume_fraction, alpha=0.5):
    """Generalized power-law dielectric mixing formula, Eq. (4.45).

    ``eps_m^alpha = eps_h^alpha + vi * (eps_i^alpha - eps_h^alpha)``

    Special cases named in Sec. 4-4.4: ``alpha=1`` is the *linear* model,
    ``alpha=1/2`` is the *refractive* model (``eps**0.5 = n`` is the
    refractive index), and ``alpha=1/3`` is the *cubic* model. This formula
    is used in Sec. 4-8.2 to build the Dobson et al. (1985) wet-soil model
    (see :func:`dobson_wet_soil_dielectric`) with ``alpha=0.65``.
    """
    host = np.asarray(eps_host, dtype=complex)
    inclusion = np.asarray(eps_inclusion, dtype=complex)
    f = np.asarray(volume_fraction, dtype=float)
    if np.any((f < 0) | (f > 1)):
        raise ValueError("volume_fraction must lie in [0, 1]")
    if alpha <= 0:
        raise ValueError("alpha must be positive")
    host_a, inclusion_a = host ** alpha, inclusion ** alpha
    return (host_a + f * (inclusion_a - host_a)) ** (1.0 / alpha)


def solid_rock_permittivity(bulk_density_g_cm3, eps_powder_real=2.0):
    """Empirical real permittivity of dry solid rock, Eq. (4.62).

    ``eps_sr' = (eps_p')^rho_b`` where ``eps_p' = 2`` is the average
    real permittivity of powdered rock (Sec. 4-7.1) and ``rho_b`` is the
    solid-rock bulk density in g/cm^3. Ulaby et al. (1990a) found ``eps_sr'``
    to be essentially frequency- and temperature-independent in the
    microwave band, so no frequency argument is required here.
    """
    rho_b = np.asarray(bulk_density_g_cm3, dtype=float)
    if np.any(rho_b < 0):
        raise ValueError("bulk_density_g_cm3 must be non-negative")
    return np.asarray(eps_powder_real, dtype=float) ** rho_b


def solid_rock_loss_factor(freq_ghz, a, b):
    """Empirical solid-rock loss factor, Eq. (4.63).

    ``eps_sr'' = a + b / f``, with rock-type-specific constants ``a`` and
    ``b`` fit from measured spectra (Ulaby et al., 1990a; Fig. 4-28). Unlike
    ``eps_sr'``, the loss factor has a definite frequency dependence
    attributed in part to bound water within the rock.
    """
    f = np.asarray(freq_ghz, dtype=float)
    if np.any(f <= 0):
        raise ValueError("freq_ghz must be positive")
    return np.asarray(a, dtype=float) + np.asarray(b, dtype=float) / f


def dry_soil_permittivity(bulk_density_g_cm3):
    """Real permittivity of dry soil, Eq. (4.64).

    ``eps_soil' = (1 + 0.44 * rho_b)^2``, essentially independent of
    frequency and temperature for dry soil (Dobson et al., 1985), with
    ``2 <~ eps_soil' <~ 4`` and ``eps_soil'' < 0.05`` typical. The loss
    factor is not modeled here because Ulaby & Long report it stays below
    0.05 regardless of density, frequency, or temperature.
    """
    rho_b = np.asarray(bulk_density_g_cm3, dtype=float)
    if np.any(rho_b < 0):
        raise ValueError("bulk_density_g_cm3 must be non-negative")
    return (1.0 + 0.44 * rho_b) ** 2


def dobson_wet_soil_dielectric(mv, freq_ghz, sand, clay, bulk_density_g_cm3=1.7,
                               temperature_c=23.0, band="wideband"):
    """Dobson et al. (1985) / Peplinski et al. (1995) semi-empirical wet-soil
    dielectric model, Eqs. (4.65)-(4.70).

    Four-component mixture (solid soil, air, bound water, free water) built
    from the generalized mixing formula (Eq. 4.45) with ``alpha=0.65``:

    ``eps_soil' = [1 + 0.66*rho_b + mv^beta1 * (eps_w')^alpha - mv]^(1/alpha)``

    ``eps_soil'' = mv^beta2 * eps_w''``

    where ``eps_w = eps_w' + j*eps_w''`` is the single-Debye water dielectric
    constant (:func:`navasar.behari.debye_water`) plus an added ionic
    conductivity loss term, Eq. (4.67b). The texture-dependent exponents are

    ``beta1 = 1.27 - 0.519*S - 0.152*C``, ``beta2 = 2.06 - 0.928*S - 0.255*C``

    (Eq. 4.68b-c), and the effective conductivity ``sigma`` [S/m] takes one
    of two empirically fit forms selected by ``band``:

    - ``"wideband"`` (1.4-18 GHz, Eq. 4.68d):
      ``sigma = -1.645 + 1.939*rho_b - 2.256*S + 1.594*C``
    - ``"lowband"`` (0.3-1.3 GHz, Peplinski et al. 1995, Eq. 4.70):
      ``sigma = 0.0467 + 0.22*rho_b - 0.411*S + 0.661*C``

    Parameters
    ----------
    mv : array_like
        Volumetric soil moisture, cm^3/cm^3.
    freq_ghz : array_like
        Frequency in GHz.
    sand, clay : float
        Sand (S) and clay (C) mass fractions in [0, 1] (silt is the
        remainder). Reported valid for 5 soil types spanning sandy loam to
        silty clay.
    bulk_density_g_cm3 : float
        Soil bulk density in g/cm^3. Ulaby & Long recommend 1.7 g/cm^3 when
        unknown.
    temperature_c : float
        Water temperature in degrees Celsius (used by the single-Debye
        water model).
    band : {"wideband", "lowband"}
        Selects the Eq. (4.68d) 1.4-18 GHz conductivity or the Peplinski
        et al. (1995) Eq. (4.70) 0.3-1.3 GHz conductivity.

    Returns
    -------
    eps : complex ndarray
        Complex relative dielectric constant of wet soil.
    """
    mv = np.asarray(mv, dtype=float)
    f = np.asarray(freq_ghz, dtype=float)
    sand, clay = float(sand), float(clay)
    rho_b = float(bulk_density_g_cm3)
    if np.any(mv < 0):
        raise ValueError("mv must be non-negative")
    if np.any(f <= 0):
        raise ValueError("freq_ghz must be positive")
    if not (0.0 <= sand <= 1.0) or not (0.0 <= clay <= 1.0) or sand + clay > 1.0:
        raise ValueError("sand and clay must be fractions in [0, 1] with sand + clay <= 1")
    if rho_b <= 0:
        raise ValueError("bulk_density_g_cm3 must be positive")

    alpha = 0.65  # Eq. (4.68a)
    beta1 = 1.27 - 0.519 * sand - 0.152 * clay  # Eq. (4.68b)
    beta2 = 2.06 - 0.928 * sand - 0.255 * clay  # Eq. (4.68c)
    if band == "wideband":
        sigma = -1.645 + 1.939 * rho_b - 2.256 * sand + 1.594 * clay  # Eq. (4.68d)
    elif band == "lowband":
        sigma = 0.0467 + 0.22 * rho_b - 0.411 * sand + 0.661 * clay  # Eq. (4.70)
    else:
        raise ValueError("band must be 'wideband' or 'lowband'")

    # Eq. (4.67a-b): single-Debye water plus added ionic conductivity loss.
    eps_w = debye_water(f, T_celsius=temperature_c)
    omega = 2.0 * np.pi * f * 1e9  # rad/s (f in GHz -> Hz)
    mv_safe = np.where(mv > 0, mv, 1.0)  # avoid divide-by-zero; masked out below
    conductivity_loss = (2.65 - rho_b) / (2.65 * mv_safe) * sigma / (_EPS0 * omega)
    eps_w_imag = eps_w.imag + np.where(mv > 0, conductivity_loss, 0.0)

    eps_soil_real = (1.0 + 0.66 * rho_b + mv ** beta1 * eps_w.real ** alpha - mv) ** (1.0 / alpha)
    eps_soil_imag = mv ** beta2 * eps_w_imag
    return eps_soil_real + 1j * eps_soil_imag
