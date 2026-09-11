# `navasar.behari` — Dielectric Models from Behari (2005)

Functions implementing dielectric, physical, radiometric, and mixing models from:
> **Behari, J. (2005).** *Microwave Dielectric Behavior of Wet Soils.* Springer Science & Business Media.

---

## Water Dielectric Models (Chapter 2)

### `debye_water(freq_ghz, T_celsius=20.0)`
Single-relaxation Debye equation for liquid water:
$$\varepsilon = \varepsilon_\infty + \frac{\varepsilon_s - \varepsilon_\infty}{1 + j\omega\tau}$$

### `cole_cole_water(freq_ghz, T_celsius=20.0, alpha=0.012)`
Cole-Cole equation accounting for relaxation time dispersion.

### `saline_water_dielectric(freq_ghz, T_celsius=20.0, salinity_ppt=35.0)`
Complex permittivity of saline water incorporating static permittivity depression and Stogryn (1971) ionic conductivity.

---

## Radiometric & Emissivity Models (Chapter 1)

### `emissivity_smooth(eps)`
Normal incidence emissivity for an ideal flat soil surface: $e_s = 1 - |(1 - \sqrt{\varepsilon}) / (1 + \sqrt{\varepsilon})|^2$.

### `emissivity_rough(e_smooth, h_roughness=0.1)`
Choudhury et al. (1982) roughness-corrected emissivity: $e_r = 1 + (e_s - 1) e^h$.

### `brightness_temperature(eps, T_soil=300.0, T_sky=5.0, T_atm=250.0, tau_atm=0.99, h_roughness=0.0)`
Radiative transfer calculation of microwave brightness temperature $T_B$.

---

## Soil Dielectric Mixing Models (Chapter 6)

### `wang_schmugge(mv, freq_ghz, sand=0.51, clay=0.13, rho_b=1.3)`
Two-regime bound-water and free-water mixing model.

### `dobson_semiemp(mv, freq_ghz, sand=0.51, clay=0.13, rho_b=1.3)`
Refractive power-law mixing model with exponent $\alpha = 0.65$.

### `rayleigh_mixing(eps_0, eps_1, f1)`
Rayleigh classical spherical inclusion model.

### `bottcher_mixing(eps_0, eps_1, f1)`
Böttcher self-consistent symmetric effective medium theory.

### `refractive_mixing(eps_list, frac_list)`
Complex Refractive Index Model (CRIM).

### `four_component_mixing(mv, freq_ghz, sand=0.51, clay=0.13, rho_b=1.3)`
Full four-component partition (solids, air, bound water, free water).
