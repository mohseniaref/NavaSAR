"""Polarimetric radar algebra and statistics based on Long & Ulaby, Chapter 5."""
import numpy as np


def fsa_to_bsa(S):
    """Convert forward-scatter alignment to backscatter alignment.

    NavaSAR uses the common convention in which the receive V basis changes
    sign on conversion: ``S_bsa = diag(1, -1) @ S_fsa``. Applying the same
    operation converts BSA back to FSA. Always record the convention with data.
    """
    matrix = np.asarray(S, dtype=complex)
    if matrix.shape[-2:] != (2, 2):
        raise ValueError("S must end in shape (2, 2)")
    transform = np.diag([1.0, -1.0])
    return np.einsum('ij,...jk->...ik', transform, matrix)


def radar_received_power(transmit_power_w, gain_linear, wavelength_m,
                         rcs_m2, range_m, system_loss=1.0):
    """Monostatic radar equation, Chapter 5 Sec. 5-4.

    Gain and system loss are linear ratios, not dB.
    """
    pt, gain, wavelength, rcs, distance, loss = np.broadcast_arrays(
        transmit_power_w, gain_linear, wavelength_m, rcs_m2, range_m, system_loss
    )
    if np.any(pt < 0) or np.any(gain < 0) or np.any(wavelength <= 0):
        raise ValueError("power/gain must be non-negative and wavelength positive")
    if np.any(rcs < 0) or np.any(distance <= 0) or np.any(loss <= 0):
        raise ValueError("RCS must be non-negative; range and loss must be positive")
    return pt * gain**2 * wavelength**2 * rcs / ((4*np.pi)**3 * distance**4 * loss)


def distributed_target_rcs(sigma0, ground_area_m2, theta_inc_deg=0.0):
    """Convert normalized backscatter to RCS using projected illuminated area."""
    sigma0 = np.asarray(sigma0, dtype=float)
    area = np.asarray(ground_area_m2, dtype=float)
    theta = np.deg2rad(np.asarray(theta_inc_deg, dtype=float))
    if np.any(sigma0 < 0) or np.any(area < 0) or np.any(np.abs(theta) >= np.pi/2):
        raise ValueError("invalid sigma0, area, or incidence angle")
    return sigma0 * area * np.cos(theta)


def scattering_matrix(shh, shv, svh, svv):
    """Construct a (..., 2, 2) scattering matrix in the H,V basis."""
    values = np.broadcast_arrays(shh, shv, svh, svv)
    result = np.empty(values[0].shape + (2, 2), dtype=complex)
    result[..., 0, 0], result[..., 0, 1] = values[0], values[1]
    result[..., 1, 0], result[..., 1, 1] = values[2], values[3]
    return result


def polarization_vector(orientation_deg=0.0, ellipticity_deg=0.0):
    """Normalized Jones vector for polarization orientation and ellipticity."""
    psi, chi = np.deg2rad([orientation_deg, ellipticity_deg])
    if abs(chi) > np.pi / 4:
        raise ValueError("ellipticity angle must lie in [-45, 45] degrees")
    return np.array([np.cos(chi) * np.cos(psi) - 1j * np.sin(chi) * np.sin(psi),
                     np.cos(chi) * np.sin(psi) + 1j * np.sin(chi) * np.cos(psi)])


def synthesize_scattering(S, transmit, receive=None):
    """Synthesize complex scattering amplitude ``e_r^H S e_t`` (Sec. 5-11)."""
    matrix = np.asarray(S, dtype=complex)
    tx = np.asarray(transmit, dtype=complex)
    rx = tx if receive is None else np.asarray(receive, dtype=complex)
    if matrix.shape[-2:] != (2, 2) or tx.shape != (2,) or rx.shape != (2,):
        raise ValueError("S must end in (2,2); polarization vectors must have length 2")
    tx = tx / np.linalg.norm(tx)
    rx = rx / np.linalg.norm(rx)
    return np.einsum('i,...ij,j->...', np.conj(rx), matrix, tx)


def pauli_vector(S, reciprocal=True):
    """Pauli target vector ``[HH+VV, HH-VV, HV+VH]/sqrt(2)``."""
    matrix = np.asarray(S, dtype=complex)
    if matrix.shape[-2:] != (2, 2):
        raise ValueError("S must end in shape (2, 2)")
    hh, hv, vh, vv = matrix[..., 0, 0], matrix[..., 0, 1], matrix[..., 1, 0], matrix[..., 1, 1]
    cross = 2.0 * hv if reciprocal else hv + vh
    return np.stack((hh + vv, hh - vv, cross), axis=-1) / np.sqrt(2.0)


def covariance_matrix(vectors, axis=0):
    """Sample covariance ``<k k^H>`` for complex scattering vectors."""
    data = np.asarray(vectors, dtype=complex)
    data = np.moveaxis(data, axis, 0)
    if data.ndim != 2:
        raise ValueError("vectors must be a 2-D array of samples by channels")
    return np.einsum('ni,nj->ij', data, np.conj(data)) / data.shape[0]


def coherency_matrix(S, sample_axis=0, reciprocal=True):
    """Pauli-basis coherency matrix ``T3=<k_p k_p^H>`` (Sec. 5-13)."""
    return covariance_matrix(pauli_vector(S, reciprocal=reciprocal), axis=sample_axis)


def polarimetric_eigendecomposition(matrix):
    """Descending eigenvalues/eigenvectors of a Hermitian covariance matrix."""
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    hermitian = (matrix + matrix.conj().T) / 2.0
    values, vectors = np.linalg.eigh(hermitian)
    order = np.argsort(values)[::-1]
    return np.maximum(values[order], 0.0), vectors[:, order]


def entropy_anisotropy(matrix):
    """Cloude-Pottier entropy and anisotropy from a 3x3 coherency matrix."""
    values, _ = polarimetric_eigendecomposition(matrix)
    if values.size != 3:
        raise ValueError("entropy/anisotropy requires a 3x3 matrix")
    total = values.sum()
    if total == 0:
        return 0.0, 0.0
    probabilities = values / total
    nonzero = probabilities > 0
    entropy = -np.sum(probabilities[nonzero] * np.log(probabilities[nonzero])) / np.log(3.0)
    denominator = values[1] + values[2]
    anisotropy = 0.0 if denominator == 0 else (values[1] - values[2]) / denominator
    return float(entropy), float(anisotropy)


def simulate_complex_speckle(covariance, n_samples=1, rng=None):
    """Draw zero-mean circular complex Gaussian scattering vectors (Secs. 5-6/5-7)."""
    covariance = np.asarray(covariance, dtype=complex)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("covariance must be square")
    values, vectors = np.linalg.eigh((covariance + covariance.conj().T) / 2)
    if values.min() < -1e-10:
        raise ValueError("covariance must be positive semidefinite")
    root = vectors @ np.diag(np.sqrt(np.maximum(values, 0.0)))
    if rng is None:
        rng = np.random.default_rng()
    noise = (rng.standard_normal((int(n_samples), values.size))
             + 1j * rng.standard_normal((int(n_samples), values.size))) / np.sqrt(2)
    return noise @ root.T


def rayleigh_amplitude(mean_power=1.0, size=None, rng=None):
    """Draw Rayleigh-fading amplitude with prescribed mean intensity."""
    if mean_power < 0:
        raise ValueError("mean_power must be non-negative")
    if rng is None:
        rng = np.random.default_rng()
    return rng.rayleigh(scale=np.sqrt(mean_power/2.0), size=size)


def multilook_intensity(mean_power=1.0, looks=1, size=None, rng=None):
    """Draw N-look intensity (Gamma law; mean P, variance P^2/N)."""
    if mean_power < 0 or int(looks) != looks or looks < 1:
        raise ValueError("mean_power must be non-negative and looks a positive integer")
    if rng is None:
        rng = np.random.default_rng()
    return rng.gamma(int(looks), float(mean_power)/int(looks), size=size)


def coherent_incoherent_power(samples, axis=0):
    """Separate coherent and incoherent power of complex samples (Sec. 5-10).

    Returns ``(|<S>|^2, <|S-<S>|^2>)``.
    """
    data = np.asarray(samples, dtype=complex)
    mean = np.mean(data, axis=axis)
    coherent = np.abs(mean)**2
    incoherent = np.mean(np.abs(data-np.expand_dims(mean, axis))**2, axis=axis)
    return coherent, incoherent


def lee_filter(intensity, window_size=5, noise_variance=None):
    """Local-statistics Lee despeckle filter (Sec. 5-9).

    ``noise_variance`` is the additive intensity-noise variance. When omitted,
    it is estimated as the median local variance. Edges are reflected.
    """
    from scipy.ndimage import uniform_filter
    image = np.asarray(intensity, dtype=float)
    if image.ndim != 2 or np.any(image < 0):
        raise ValueError("intensity must be a non-negative 2-D array")
    if int(window_size) != window_size or window_size < 3 or window_size % 2 == 0:
        raise ValueError("window_size must be an odd integer >= 3")
    local_mean = uniform_filter(image, int(window_size), mode="reflect")
    local_second = uniform_filter(image**2, int(window_size), mode="reflect")
    local_variance = np.maximum(local_second-local_mean**2, 0.0)
    noise = np.median(local_variance) if noise_variance is None else float(noise_variance)
    if noise < 0:
        raise ValueError("noise_variance must be non-negative")
    weight = np.maximum(local_variance-noise, 0.0) / np.maximum(local_variance, 1e-15)
    return local_mean + weight*(image-local_mean)


def polarimetric_parameters(matrix):
    """Return span, eigenvalue probabilities, entropy and anisotropy.

    Accepts a 3x3 covariance or coherency matrix. These compact parameters
    describe total power and randomness; they are not a land-cover classifier.
    """
    values, _ = polarimetric_eigendecomposition(matrix)
    if values.size != 3:
        raise ValueError("polarimetric parameters require a 3x3 matrix")
    span = float(values.sum())
    probabilities = np.zeros(3) if span == 0 else values/span
    entropy, anisotropy = entropy_anisotropy(matrix)
    return {"span": span, "probabilities": probabilities,
            "entropy": entropy, "anisotropy": anisotropy}
