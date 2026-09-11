import numpy as np
import pytest

from navasar.materials import (
    brine_salinity, canopy_effective_permittivity, depolarization_factors,
    dilute_ellipsoid_mixing, dry_snow_permittivity, dry_snow_tvb,
    polder_van_santen, powdered_rock_permittivity, pure_ice_permittivity,
    refractive_mixture, sea_ice_permittivity, tvb_mixing,
    vegetation_permittivity, wet_snow_permittivity,
)


def test_sphere_depolarization():
    np.testing.assert_allclose(depolarization_factors(1, 1, 1), [1/3] * 3, atol=1e-8)


def test_ellipsoid_endpoints():
    assert dilute_ellipsoid_mixing(2, 8, 0) == pytest.approx(2)
    assert dilute_ellipsoid_mixing(2, 8, 0.1).real > 2


def test_polder_endpoints_and_symmetry():
    assert polder_van_santen([2, 8], [1, 0]) == pytest.approx(2)
    a = polder_van_santen([2, 8], [0.6, 0.4])
    b = polder_van_santen([8, 2], [0.4, 0.6])
    assert a == pytest.approx(b)


def test_dry_snow_and_refractive_limits():
    assert dry_snow_permittivity(0) == pytest.approx(1)
    assert dry_snow_permittivity(0.4).real > 1
    assert refractive_mixture([3 + 0.1j, 80 + 10j], [1, 0]) == pytest.approx(3 + 0.1j)


def test_ice_and_cryosphere_models_are_passive():
    ice = pure_ice_permittivity(np.array([1.0, 10.0]), -10)
    assert np.all(ice.real > 3) and np.all(ice.imag > 0)
    assert brine_salinity(-5) > 35
    sea_ice = sea_ice_permittivity(5.0, -10, 5.0)
    assert sea_ice.real > ice[0].real and sea_ice.imag > 0
    assert dry_snow_tvb(5, 0).real == pytest.approx(1)
    assert wet_snow_permittivity(5, 0.25, 0.05).imag > 0


def test_tvb_endpoints_and_material_models():
    assert tvb_mixing(2, 8, 0) == pytest.approx(2)
    assert tvb_mixing(2, 8, 1) == pytest.approx(8)
    assert powdered_rock_permittivity(6, 0).real == pytest.approx(1)
    plant = vegetation_permittivity(5, 0.7)
    canopy = canopy_effective_permittivity(plant, 0.02)
    assert plant.real > canopy.real > 1
