import numpy as np
import pytest

from navasar.materials import (
    brine_salinity, canopy_effective_permittivity, depolarization_factors,
    dilute_ellipsoid_mixing, dobson_wet_soil_dielectric, dry_snow_permittivity,
    dry_snow_tvb, dry_soil_permittivity, generalized_dielectric_mixing,
    polder_van_santen, powdered_rock_permittivity, pure_ice_permittivity,
    refractive_mixture, sea_ice_permittivity, solid_rock_loss_factor,
    solid_rock_permittivity, tvb_mixing, vegetation_permittivity,
    wet_snow_permittivity,
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


def test_generalized_mixing_endpoints_and_named_cases():
    # Endpoints reduce to the host/inclusion permittivity for any alpha.
    assert generalized_dielectric_mixing(2, 8, 0, alpha=0.65) == pytest.approx(2)
    assert generalized_dielectric_mixing(2, 8, 1, alpha=0.65) == pytest.approx(8)
    # alpha=1 is the linear (arithmetic) model.
    assert generalized_dielectric_mixing(2, 8, 0.3, alpha=1.0) == pytest.approx(2 + 0.3 * 6)
    # alpha=0.5 is the refractive-index model: sqrt(eps_m) is linear in vi.
    linear_index = (1 - 0.3) * np.sqrt(2) + 0.3 * np.sqrt(8)
    assert generalized_dielectric_mixing(2, 8, 0.3, alpha=0.5) == pytest.approx(linear_index ** 2)
    with pytest.raises(ValueError):
        generalized_dielectric_mixing(2, 8, 0.5, alpha=0)


def test_solid_rock_and_dry_soil_density_models():
    # Eq. (4.62): eps_sr' = eps_p'^rho_b, eps_p' = 2 by default.
    assert solid_rock_permittivity(0) == pytest.approx(1.0)
    assert solid_rock_permittivity(1) == pytest.approx(2.0)
    assert solid_rock_permittivity(2.5) > solid_rock_permittivity(1.0)
    with pytest.raises(ValueError):
        solid_rock_permittivity(-1)

    # Eq. (4.63): eps_sr'' = a + b/f, decreasing with frequency for b > 0.
    loss = solid_rock_loss_factor(np.array([1.0, 10.0]), a=0.01, b=0.05)
    assert loss[0] > loss[1] > 0
    with pytest.raises(ValueError):
        solid_rock_loss_factor(0, a=0.01, b=0.05)

    # Eq. (4.64): eps_soil' = (1 + 0.44*rho_b)^2, typical range 2-4.
    assert dry_soil_permittivity(0) == pytest.approx(1.0)
    dry = dry_soil_permittivity(1.7)
    assert 2.0 < dry < 4.0
    with pytest.raises(ValueError):
        dry_soil_permittivity(-1)


def test_dobson_wet_soil_dielectric_matches_book_behavior():
    mv = np.array([0.0, 0.1, 0.2, 0.3])
    eps = dobson_wet_soil_dielectric(mv, freq_ghz=1.4, sand=0.5, clay=0.15,
                                     bulk_density_g_cm3=1.7, band="wideband")
    # Both parts should increase monotonically with moisture (Fig. 4-33/4-34).
    assert np.all(np.diff(eps.real) > 0)
    assert np.all(np.diff(eps.imag[1:]) > 0)  # imag[0] is exactly 0 at mv=0
    assert eps.imag[0] == pytest.approx(0.0)

    # The 0.3-1.3 GHz (Peplinski) conductivity differs from the 1.4-18 GHz
    # (Dobson) conductivity, so the loss factor differs for mv > 0 while the
    # real part (independent of sigma) matches exactly.
    eps_low = dobson_wet_soil_dielectric(mv, freq_ghz=1.0, sand=0.5, clay=0.15,
                                         bulk_density_g_cm3=1.7, band="lowband")
    eps_wide = dobson_wet_soil_dielectric(mv, freq_ghz=1.0, sand=0.5, clay=0.15,
                                          bulk_density_g_cm3=1.7, band="wideband")
    np.testing.assert_allclose(eps_low.real, eps_wide.real)
    assert not np.allclose(eps_low.imag[1:], eps_wide.imag[1:])

    with pytest.raises(ValueError):
        dobson_wet_soil_dielectric(0.2, 1.4, sand=0.7, clay=0.5)  # sand+clay > 1
    with pytest.raises(ValueError):
        dobson_wet_soil_dielectric(0.2, 1.4, sand=0.5, clay=0.15, band="bad")
