# Land backscatter and SAR image simulation

`navasar.land` implements the reusable quantitative parts of Chapter 21 in
Ulaby, Moore & Fung, *Microwave Remote Sensing: Active and Passive*, volume
III. The implementation is based on the printed equations.

This module works with calibrated, detected power (`sigma0`, linear units).
It is intentionally separate from `navasar.slc`, which produces complex,
phase-preserving soil echoes for interferometry.

## Implemented model chain

For a vegetation canopy above a rough ground surface, the first-order model is

$$
\sigma^0_\mathrm{can} = \sigma^0_v + T^2\sigma^0_s + \sigma^0_{sv},
\qquad
T^2 = \exp(-2\tau\sec\theta),
$$

where $\tau$ is one-way optical depth. With a Rayleigh phase function and
single scattering (co-polarization),

$$
\sigma^0_v = \frac{3}{4}\omega
  [1-\exp(-2\tau\sec\theta)]\cos\theta.
$$

These are Eqs. (21.50)–(21.52). Use `canopy_transmissivity`,
`vegetation_volume_backscatter`, and `canopy_backscatter`.

`periodic_surface_backscatter` numerically evaluates the Eq. (21.37) facet
integral for sinusoidal agricultural rows. The small-scale surface law is a
callable, so a validated scattering model can be supplied independently.

For a detected $N$-look intensity image, `simulate_multilook_intensity` uses
the normalized Gamma fading law (mean 1, variance $1/N$). With independent,
unit-mean texture variance $\sigma_T^2$, Eq. (21.89a) gives

$$
\frac{\mathrm{var}(I)}{E[I]^2}
=\sigma_T^2+\frac{1}{N}+\frac{\sigma_T^2}{N}.
$$

`intensity_moments` evaluates that relation and
`estimate_texture_variance` implements its inversion, Eq. (21.89b).

## Example

```python
import numpy as np
from navasar.land import canopy_backscatter, simulate_multilook_intensity

bare_soil = np.full((256, 256), 0.08)  # linear sigma0 (about -11 dB)
mean_land = canopy_backscatter(
    bare_soil, optical_depth=0.35, albedo=0.12, theta_inc_deg=30
)
image = simulate_multilook_intensity(
    mean_land, looks=4, rng=np.random.default_rng(7)
)
```

## Scope and limitations

- The historical land-application examples and empirical observations in the
  chapter are documentation, not universal forward models.
- The optional surface–volume interaction term is accepted explicitly. The
  chapter's fitted interaction formulas have narrow parameter ranges and are
  therefore not silently applied outside those experiments.
- Cross-polarized volume backscatter is zero only under the first-order
  Rayleigh assumptions used by Eq. (21.52); real vegetation can have a strong
  cross-polarized return.
