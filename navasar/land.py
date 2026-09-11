"""Land backscatter and detected-image simulation.

The equations in this module follow Chapter 21 of Ulaby, Moore & Fung,
*Microwave Remote Sensing: Active and Passive*, vol. III (1986).  Inputs and
outputs called ``sigma0`` are dimensionless linear power ratios, not decibels.

This module deliberately separates detected land imagery from
``navasar.slc``: the latter simulates complex, phase-preserving soil echoes,
whereas Chapter 21 primarily describes calibrated backscatter and image
intensity statistics.
"""

import numpy as np


def db_to_linear(value_db):
    """Convert a power quantity in dB to a dimensionless linear ratio."""
    return np.power(10.0, np.asarray(value_db, dtype=float) / 10.0)


def linear_to_db(value, floor_db=None):
    """Convert a non-negative power ratio to dB.

    ``floor_db`` can be used to replace ``-inf`` for zero-valued pixels.
    """
    value = np.asarray(value, dtype=float)
    if np.any(value < 0):
        raise ValueError("power values must be non-negative")
    with np.errstate(divide="ignore"):
        result = 10.0 * np.log10(value)
    if floor_db is not None:
        result = np.maximum(result, float(floor_db))
    return result


def canopy_transmissivity(optical_depth, theta_inc_deg=30.0):
    """Two-way vegetation transmissivity, Chapter 21 Eq. (21.50).

    Returns ``exp(-2 * optical_depth * sec(theta))``.
    """
    tau = np.asarray(optical_depth, dtype=float)
    theta = np.deg2rad(np.asarray(theta_inc_deg, dtype=float))
    if np.any(tau < 0):
        raise ValueError("optical_depth must be non-negative")
    if np.any(np.abs(theta) >= np.pi / 2):
        raise ValueError("incidence angle must lie between -90 and 90 degrees")
    return np.exp(-2.0 * tau / np.cos(theta))


def vegetation_volume_backscatter(optical_depth, albedo,
                                  theta_inc_deg=30.0):
    """Rayleigh, single-scattering canopy term, Eq. (21.52).

    This is the co-polarized volume contribution.  Under the assumptions of
    Eq. (21.52), the cross-polarized contribution is zero.
    """
    tau = np.asarray(optical_depth, dtype=float)
    omega = np.asarray(albedo, dtype=float)
    if np.any(tau < 0):
        raise ValueError("optical_depth must be non-negative")
    if np.any((omega < 0) | (omega > 1)):
        raise ValueError("albedo must be between 0 and 1")
    theta = np.deg2rad(np.asarray(theta_inc_deg, dtype=float))
    if np.any(np.abs(theta) >= np.pi / 2):
        raise ValueError("incidence angle must lie between -90 and 90 degrees")
    return 0.75 * omega * (1.0 - np.exp(-2.0 * tau / np.cos(theta))) * np.cos(theta)


def canopy_backscatter(surface_sigma0, optical_depth, albedo,
                       theta_inc_deg=30.0, interaction_sigma0=0.0):
    """First-order canopy-over-ground backscatter, Eqs. (21.50)-(21.52).

    ``sigma0_can = sigma0_volume + T2*sigma0_surface + sigma0_interaction``.
    The optional interaction term is explicit because the textbook's fitted
    interaction models have narrow validity ranges and are not universal.
    """
    surface = np.asarray(surface_sigma0, dtype=float)
    interaction = np.asarray(interaction_sigma0, dtype=float)
    if np.any(surface < 0) or np.any(interaction < 0):
        raise ValueError("backscatter terms must be non-negative")
    volume = vegetation_volume_backscatter(
        optical_depth, albedo, theta_inc_deg=theta_inc_deg
    )
    return (volume
            + canopy_transmissivity(optical_depth, theta_inc_deg) * surface
            + interaction)


def periodic_surface_backscatter(local_sigma0, theta_inc_deg=30.0,
                                 look_azimuth_deg=90.0, amplitude=0.05,
                                 wavelength=0.60, n_samples=2048):
    """Integrate local backscatter over a sinusoidal row surface (Eq. 21.37).

    Parameters
    ----------
    local_sigma0 : callable
        Function accepting local incidence angle in degrees and returning
        linear local backscatter. It supplies the small-scale roughness model.
    look_azimuth_deg : float
        Look azimuth measured from the x axis; rows vary along y. Thus 90 deg
        is perpendicular to rows and 0 deg is parallel.
    amplitude, wavelength : float
        Sinusoid amplitude and spatial period in the same length unit.

    Notes
    -----
    The returned mean contains the surface-area factor ``sec(alpha)``. Samples
    whose local facet faces away from the radar are excluded.
    """
    if amplitude < 0 or wavelength <= 0:
        raise ValueError("amplitude must be >= 0 and wavelength must be > 0")
    if int(n_samples) < 16:
        raise ValueError("n_samples must be at least 16")
    theta = np.deg2rad(float(theta_inc_deg))
    phi = np.deg2rad(float(look_azimuth_deg))
    if not 0 <= theta < np.pi / 2:
        raise ValueError("theta_inc_deg must be in [0, 90)")

    y = (np.arange(int(n_samples)) + 0.5) * wavelength / int(n_samples)
    slope = (2.0 * np.pi * amplitude / wavelength) * np.cos(2.0 * np.pi * y / wavelength)
    sec_alpha = np.sqrt(1.0 + slope ** 2)
    cos_local = (np.cos(theta) - slope * np.sin(theta) * np.sin(phi)) / sec_alpha
    visible = cos_local > 0
    local_angle = np.rad2deg(np.arccos(np.clip(cos_local[visible], 0.0, 1.0)))
    sigma = np.asarray(local_sigma0(local_angle), dtype=float)
    if sigma.shape not in ((), local_angle.shape):
        raise ValueError("local_sigma0 must return a scalar or one value per angle")
    if np.any(sigma < 0):
        raise ValueError("local_sigma0 returned negative backscatter")
    contributions = np.zeros_like(cos_local)
    contributions[visible] = sigma * sec_alpha[visible]
    return contributions.mean()


def intensity_moments(mean_sigma0, looks, texture_variance=0.0):
    """Mean and variance of the multiplicative intensity model (Eq. 21.89a).

    ``texture_variance`` is the variance of a unit-mean texture variable.
    """
    mean = np.asarray(mean_sigma0, dtype=float)
    n = np.asarray(looks, dtype=float)
    tex_var = np.asarray(texture_variance, dtype=float)
    if np.any(mean < 0) or np.any(n <= 0) or np.any(tex_var < 0):
        raise ValueError("mean and texture variance must be non-negative; looks must be positive")
    normalized_variance = tex_var + 1.0 / n + tex_var / n
    return mean, mean ** 2 * normalized_variance


def estimate_texture_variance(intensity, looks):
    """Estimate unit-mean texture variance from an image, Eq. (21.89b)."""
    image = np.asarray(intensity, dtype=float)
    if image.size < 2 or np.any(image < 0):
        raise ValueError("intensity must contain at least two non-negative samples")
    if looks <= 0:
        raise ValueError("looks must be positive")
    mean = image.mean()
    if mean == 0:
        return 0.0
    observed = image.var(ddof=1) / mean ** 2
    return max(0.0, (observed - 1.0 / looks) / (1.0 + 1.0 / looks))


def simulate_multilook_intensity(mean_sigma0, looks=1, texture=None, rng=None):
    """Simulate an N-look detected SAR image using Chapter 21 statistics.

    Rayleigh fading after incoherent averaging is Gamma distributed with
    shape ``looks`` and scale ``1/looks`` (unit mean, variance ``1/looks``).
    ``texture``, when supplied, is a non-negative unit-mean multiplier (or a
    broadcastable array of multipliers).
    """
    mean = np.asarray(mean_sigma0, dtype=float)
    if np.any(mean < 0):
        raise ValueError("mean_sigma0 must be non-negative")
    if int(looks) != looks or looks < 1:
        raise ValueError("looks must be a positive integer")
    if rng is None:
        rng = np.random.default_rng()
    if texture is None:
        texture = 1.0
    texture = np.asarray(texture, dtype=float)
    if np.any(texture < 0):
        raise ValueError("texture must be non-negative")
    fading = rng.gamma(shape=int(looks), scale=1.0 / int(looks), size=mean.shape)
    return mean * texture * fading

