import numpy as np
import pytest

from navasar.land import (
    canopy_backscatter,
    canopy_transmissivity,
    db_to_linear,
    estimate_texture_variance,
    intensity_moments,
    linear_to_db,
    periodic_surface_backscatter,
    simulate_multilook_intensity,
    vegetation_volume_backscatter,
)


def test_db_round_trip():
    values = np.array([-30.0, -10.0, 0.0, 7.5])
    np.testing.assert_allclose(linear_to_db(db_to_linear(values)), values)


def test_canopy_limits_and_composition():
    surface = 0.2
    assert canopy_transmissivity(0.0, 30.0) == pytest.approx(1.0)
    assert vegetation_volume_backscatter(0.0, 0.2, 30.0) == pytest.approx(0.0)
    assert canopy_backscatter(surface, 0.0, 0.2, 30.0) == pytest.approx(surface)
    assert canopy_transmissivity(0.8, 30.0) < canopy_transmissivity(0.2, 30.0)


def test_flat_periodic_surface_reduces_to_local_law():
    law = lambda theta: np.full_like(theta, 0.125)
    result = periodic_surface_backscatter(law, amplitude=0.0)
    assert result == pytest.approx(0.125)


def test_rows_create_look_direction_modulation():
    law = lambda theta: np.exp(-np.square(theta / 18.0))
    parallel = periodic_surface_backscatter(
        law, look_azimuth_deg=0.0, amplitude=0.055, wavelength=0.606
    )
    perpendicular = periodic_surface_backscatter(
        law, look_azimuth_deg=90.0, amplitude=0.055, wavelength=0.606
    )
    assert parallel != pytest.approx(perpendicular)


def test_multilook_statistics_match_theory():
    rng = np.random.default_rng(1234)
    mean = np.full(300_000, 0.2)
    image = simulate_multilook_intensity(mean, looks=4, rng=rng)
    expected_mean, expected_var = intensity_moments(0.2, 4)
    assert image.mean() == pytest.approx(expected_mean, rel=0.01)
    assert image.var() == pytest.approx(expected_var, rel=0.02)
    assert estimate_texture_variance(image, looks=4) == pytest.approx(0.0, abs=0.01)


def test_invalid_inputs():
    with pytest.raises(ValueError):
        canopy_transmissivity(-0.1)
    with pytest.raises(ValueError):
        simulate_multilook_intensity(np.ones(2), looks=0)
