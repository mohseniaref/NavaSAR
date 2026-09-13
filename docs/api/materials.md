# Material dielectric models

`navasar.materials` collects reusable Chapter 4 dielectric building blocks.
It reuses the existing water implementations and adds ice, snow, sea ice,
rock, soil, vegetation, and general effective-medium relations.

Model groups include water and ice, dry and wet snow, sea ice, solid and
powdered rock, dry and wet soil, vegetation material and canopy-scale
mixtures. General mixture tools include ellipsoid depolarization factors,
dilute mixing, Polder–van Santen, TVB, the generalized power-law formula,
and refractive-index formulations.

Complex values use the NavaSAR loss convention $\epsilon=\epsilon'+j\epsilon''$.
Mixture fractions are volume fractions and must sum to one where applicable.
The compact sea-ice and wet-snow functions remain effective-medium models and
do not replace ice-type-specific field calibration.

## Generalized dielectric mixing (Sec. 4-4.4)

`generalized_dielectric_mixing(eps_host, eps_inclusion, volume_fraction, alpha)`
implements the power-law formula of Eq. (4.45),

$$
\varepsilon_m^\alpha = \varepsilon_h^\alpha + v_i\,(\varepsilon_i^\alpha - \varepsilon_h^\alpha),
$$

which unifies three named special cases: $\alpha=1$ is the *linear* model,
$\alpha=1/2$ is the *refractive* model (since $\varepsilon^{1/2}=n$ is the
refractive index — this is the two-component version of
`refractive_mixture`), and $\alpha=1/3$ is the *cubic* model. Section 4-8.2
uses this same formula with $\alpha=0.65$ to build the Dobson wet-soil model
below.

## Solid rock (Sec. 4-7.2)

- `solid_rock_permittivity(bulk_density_g_cm3, eps_powder_real=2.0)` —
  Eq. (4.62): $\varepsilon_{sr}'=(\varepsilon_p')^{\rho_b}$. Ulaby et al.
  (1990a) found $\varepsilon_{sr}'$ essentially frequency- and
  temperature-independent across the microwave band.
- `solid_rock_loss_factor(freq_ghz, a, b)` — Eq. (4.63):
  $\varepsilon_{sr}''=a+b/f$, with rock-type-specific constants $a,b$.
  Unlike $\varepsilon_{sr}'$, the loss factor varies clearly with frequency.

## Soil (Sec. 4-8)

- `dry_soil_permittivity(bulk_density_g_cm3)` — Eq. (4.64):
  $\varepsilon_{soil}'=(1+0.44\rho_b)^2$ for dry soil, essentially
  frequency- and temperature-independent, with $\varepsilon_{soil}''<0.05$.
- `dobson_wet_soil_dielectric(mv, freq_ghz, sand, clay, bulk_density_g_cm3=1.7,
  temperature_c=23.0, band="wideband")` — the full Dobson et al. (1985)
  four-component semi-empirical model, Eqs. (4.65)-(4.69), together with the
  Peplinski et al. (1995) low-band conductivity correction, Eq. (4.70). This
  is a more literal, texture-resolved implementation than the compact
  `navasar.behari.dobson_semiemp`: it derives $\beta_1,\beta_2$ from sand/clay
  mass fractions (Eq. 4.68b-c) and lets `band` select between the 1.4-18 GHz
  conductivity (Eq. 4.68d, `"wideband"`) and the 0.3-1.3 GHz conductivity
  (Eq. 4.70, `"lowband"`) used at L-band and below.

```python
import numpy as np
from navasar.materials import dobson_wet_soil_dielectric

mv = np.linspace(0.0, 0.4, 5)
eps_lband = dobson_wet_soil_dielectric(
    mv, freq_ghz=1.4, sand=0.51, clay=0.13, band="lowband"
)
```
