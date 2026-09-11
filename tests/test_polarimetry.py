import numpy as np
import pytest

from navasar.polarimetry import (
    coherent_incoherent_power, covariance_matrix, distributed_target_rcs,
    entropy_anisotropy, fsa_to_bsa, lee_filter, multilook_intensity,
    pauli_vector, polarimetric_parameters, polarization_vector,
    radar_received_power, rayleigh_amplitude, scattering_matrix,
    simulate_complex_speckle, synthesize_scattering,
)


def test_linear_polarization_synthesis():
    S = scattering_matrix(2, 0.2, 0.2, 1)
    h = polarization_vector(0, 0)
    v = polarization_vector(90, 0)
    assert synthesize_scattering(S, h) == pytest.approx(2)
    assert synthesize_scattering(S, v) == pytest.approx(1)


def test_pauli_vector_and_covariance_are_hermitian():
    S = scattering_matrix([1, 2], [0, 0.2j], [0, 0.2j], [1, 0.5])
    k = pauli_vector(S)
    C = covariance_matrix(k)
    np.testing.assert_allclose(C, C.conj().T)
    assert np.linalg.eigvalsh(C).min() > -1e-12


def test_entropy_bounds():
    H, A = entropy_anisotropy(np.diag([3.0, 2.0, 1.0]))
    assert 0 <= H <= 1
    assert 0 <= A <= 1


def test_complex_speckle_recovers_covariance():
    target = np.array([[2, 0.3 + 0.2j], [0.3 - 0.2j, 1]])
    samples = simulate_complex_speckle(target, 250_000, np.random.default_rng(5))
    measured = covariance_matrix(samples)
    np.testing.assert_allclose(measured, target, rtol=0.015, atol=0.015)


def test_convention_and_radar_equation():
    S = scattering_matrix(1, 2, 3, 4)
    np.testing.assert_allclose(fsa_to_bsa(fsa_to_bsa(S)), S)
    rcs = distributed_target_rcs(0.1, 100, 60)
    assert rcs == pytest.approx(5)
    assert radar_received_power(1000, 100, 0.1, rcs, 1000) > 0


def test_fading_coherent_split_and_filter():
    rng = np.random.default_rng(8)
    amp = rayleigh_amplitude(2, 200_000, rng)
    assert np.mean(amp**2) == pytest.approx(2, rel=.01)
    looks = multilook_intensity(2, 4, 200_000, rng)
    assert looks.var() == pytest.approx(1, rel=.02)
    coherent, incoherent = coherent_incoherent_power(np.ones(100)+.1j)
    assert coherent == pytest.approx(1.01) and incoherent == pytest.approx(0)
    filtered = lee_filter(rng.gamma(1, 1, (50, 50)))
    assert filtered.var() < 1.0
    assert polarimetric_parameters(np.eye(3))["entropy"] == pytest.approx(1)
