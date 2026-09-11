# `navasar.coherence` — Interferometric Coherence and Phase Triplets

**Source:** `navasar/coherence.py`

This module answers: **given two soil moisture values, what is the
interferometric phase and coherence between the two SAR images?**

---

## Background: what is interferometric coherence?

When you take two SAR images of the same area at different times, you
can multiply them together pixel by pixel to form an **interferogram**.
The phase of each pixel in the interferogram is the phase *difference*
between the two images.

If nothing changed between the two acquisitions, the phase difference
is zero everywhere and the images are perfectly **coherent** (|γ| = 1).
If the soil moisture changed, the wave penetrates to a different depth
and accumulates a different phase — the images are partially decorrelated
(|γ| < 1) and there is a systematic phase bias (∠γ ≠ 0).

---

## `interferometric_coherence`

```python
interferometric_coherence(mv1, mv2, freq_ghz=1.4, theta_inc_deg=45.0,
                           sand=0.51, clay=0.13)
```

### What it does

Computes the complex coherence γ between a master image (moisture mv1)
and a slave image (moisture mv2), using the analytical model of
De Zan et al. (2014), Eq. 12–14.

### The physics behind the formula

**Step 1 — model a SAR pixel as a depth integral.**

Each pixel is the sum of contributions from scatterers at all depths z:

```
p = ∫₀^∞  ξ(z) · exp(−2j·kz'·z) dz
```

- `ξ(z)` = random scatterer at depth z (we don't know it, but its
  statistics are known)
- `exp(−2j·kz'·z)` = two-way propagation phasor to depth z
- The factor 2 accounts for the wave going down AND coming back up

**Step 2 — compute the expected interferogram.**

The expected value of p₁·p₂* (master × slave conjugate) is:

```
I(ε₁, ε₂) = E[p₁·p₂*]
           = ∫₀^∞ f(z) · exp(−2j·kz1'·z) · exp(+2j·kz2'*·z) dz
```

For a **uniform scatterer density** f(z) = constant, this integral
has a closed-form solution:

```
I(ε₁, ε₂) = 1 / (2j·kz1' − 2j·kz2'*)
```

**Why does this converge?** The integrand is exp(−(2j·kz1' − 2j·kz2'*)·z).
Since Im(kz') < 0, the real part of the exponent is negative, so the
integral converges to 1/(2j·kz1' − 2j·kz2'*).

**Step 3 — normalise to get coherence.**

The coherence normalises by the power of each image:

```
γ = I(ε₁, ε₂) / sqrt(I(ε₁,ε₁) · I(ε₂,ε₂))
```

where:
```
I(ε₁,ε₁) = 1 / (2j·kz1' − 2j·kz1'*)
           = 1 / (−4·Im(kz1'))      ← purely real, positive
I(ε₂,ε₂) = 1 / (−4·Im(kz2'))      ← purely real, positive
```

### What the output means

```
γ = |γ| · exp(j·φ)
     ↑          ↑
  coherence   interferometric
  magnitude   phase bias
```

**Coherence magnitude |γ|:**
- = 1.0 when mv1 = mv2 (no change → perfect coherence)
- < 1.0 when mv1 ≠ mv2 (moisture changed → decorrelation)
- Decreases as |mv1 − mv2| increases

**Interferometric phase ∠γ:**
- = 0 when mv1 = mv2
- ≠ 0 when moisture changed
- Caused by the change in Re(kz') — the wave accumulates a different
  phase at each depth when moisture changes
- Typical magnitude: 10–50° for realistic moisture changes at L-band

**Important:** only a change in Re(kz') shifts the phase. A change in
Im(kz') only changes the penetration depth and therefore |γ|, not ∠γ.
This means a soil that simply dries out uniformly (only penetration
depth changes) produces coherence loss but no phase bias.

### Numerical walkthrough

Master: mv1 = 0.20, Slave: mv2 = 0.25, L-band, 45° incidence.

```
kz1' ≈ −94.5 − 11.2j  rad/m   (at mv=0.20)
kz2' ≈ −97.1 − 12.6j  rad/m   (at mv=0.25)

I(ε₁,ε₂) = 1 / (2j·(−94.5−11.2j) − 2j·(−97.1+12.6j))
          = 1 / (2j·(−94.5−11.2j) − 2j·(−97.1+12.6j))
          = 1 / (−189j + 22.4 + 194.2j − 25.2)
          = 1 / (−2.8 + 5.2j)

I(ε₁,ε₁) = 1/(4×11.2) = 0.02232
I(ε₂,ε₂) = 1/(4×12.6) = 0.01984

|γ| = |I(ε₁,ε₂)| / sqrt(I(ε₁,ε₁)·I(ε₂,ε₂))
    ≈ 0.167 / sqrt(0.02232 × 0.01984) ≈ 0.79

∠γ ≈ −7°   (phase bias due to moisture change)
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv1` | array_like | — | Master moisture. Can be a scalar, 1-D array, or 2-D grid. |
| `mv2` | array_like | — | Slave moisture. Must broadcast with mv1. |
| `freq_ghz` | float | `1.4` | Radar frequency. L-band=1.4, C-band=5.0. |
| `theta_inc_deg` | float | `45.0` | Incidence angle in degrees. |
| `sand` | float | `0.51` | Sand fraction [0,1]. |
| `clay` | float | `0.13` | Clay fraction [0,1]. |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `gamma` | complex ndarray | Complex coherence. Use `abs(gamma)` for magnitude, `np.rad2deg(np.angle(gamma))` for phase in degrees. |

### Example — reproduce De Zan (2014) Fig. 3 and Fig. 4

```python
import numpy as np
import matplotlib.pyplot as plt
from navasar.coherence import interferometric_coherence

mv = np.linspace(0.0, 0.5, 100)
MV1, MV2 = np.meshgrid(mv, mv)   # master on x-axis, slave on y-axis
gamma = interferometric_coherence(MV1, MV2, freq_ghz=1.4)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
im0 = axes[0].contourf(mv, mv, np.abs(gamma), levels=20, cmap='viridis')
plt.colorbar(im0, ax=axes[0], label='|γ|')
axes[0].set_title('Fig. 3 — Coherence magnitude')

im1 = axes[1].contourf(mv, mv, np.rad2deg(np.angle(gamma)), levels=20, cmap='RdBu_r')
plt.colorbar(im1, ax=axes[1], label='Phase [deg]')
axes[1].set_title('Fig. 4 — Interferometric phase')

for ax in axes:
    ax.set_xlabel('Master moisture mv1')
    ax.set_ylabel('Slave moisture mv2')
plt.tight_layout()
```

---

## `phase_triplet`

```python
phase_triplet(mv1, mv2, mv3, freq_ghz=1.4, theta_inc_deg=45.0,
              sand=0.51, clay=0.13)
```

### What it does

Computes the **phase triplet** (also called closure phase) for three
acquisitions. This is the key observable that reveals soil moisture
changes without needing absolute phase calibration.

### Why triplets are useful

The interferometric phase ∠γ contains the moisture signal, but also:
- Orbital errors (constant offset per image)
- Atmospheric delay (varies per image)
- Deformation (varies per image)

These nuisance terms are **linear** — they add a constant to each image's
phase. The triplet cancels them automatically.

### The formula

Given three acquisitions with moistures mv1, mv2, mv3:

```
φ₁₂ = ∠γ(mv1, mv2)   ← phase of interferogram 1→2
φ₂₃ = ∠γ(mv2, mv3)   ← phase of interferogram 2→3
φ₁₃ = ∠γ(mv1, mv3)   ← phase of interferogram 1→3

triplet = φ₁₂ + φ₂₃ − φ₁₃
```

Equivalently, using the triple product of coherences:

```
triplet = ∠(γ₁₂ · γ₂₃ · γ₁₃*)
```

The `*` denotes complex conjugate. This form is numerically more stable.

### Why is the triplet non-zero for moisture?

For a **linear** signal (deformation d per image):
```
φ₁₂ = d₂ − d₁
φ₂₃ = d₃ − d₂
φ₁₃ = d₃ − d₁
triplet = (d₂−d₁) + (d₃−d₂) − (d₃−d₁) = 0  ✓ always zero
```

For **soil moisture**, the phase is a non-linear function of mv.
The path 1→2→3 is not the same as 1→3 directly, because the
intermediate moisture state mv2 affects the scattering differently.

**Symmetric case (mv1 = mv3):** triplet = 0 always, regardless of mv2.
This is because the path 1→2→3 is the reverse of 3→2→1, and the
non-linearity cancels.

**Asymmetric case (mv1 ≠ mv3):** triplet ≠ 0. The larger the difference
between mv1 and mv3, and the more extreme mv2 is, the larger the triplet.

### Typical magnitudes

| Scenario | Triplet magnitude |
|----------|------------------|
| mv1=0.10, mv2=0.25, mv3=0.15, L-band | ≈ 22° |
| mv1=0.05, mv2=0.20, mv3=0.40, L-band | ≈ 40° |
| mv1=mv3=0.20, any mv2, L-band | = 0° (symmetric) |
| Same scenario at C-band | smaller (less penetration) |

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv1` | array_like | — | Moisture at acquisition 1 (master). |
| `mv2` | array_like | — | Moisture at acquisition 2 (middle). This is the one that "breaks" the linearity. |
| `mv3` | array_like | — | Moisture at acquisition 3 (slave). |
| `freq_ghz` | float | `1.4` | Radar frequency. |
| `theta_inc_deg` | float | `45.0` | Incidence angle. |
| `sand` | float | `0.51` | Sand fraction. |
| `clay` | float | `0.13` | Clay fraction. |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `triplet` | ndarray | Phase mismatch in **degrees**. Positive or negative depending on the direction of moisture change. |

### Example

```python
import numpy as np
from navasar.coherence import phase_triplet

# AGRISAR 2006 moisture series (May-11 to Jul-05)
mv_series = np.array([0.28, 0.22, 0.18, 0.15, 0.12, 0.10])

from itertools import combinations
for i, j, k in combinations(range(6), 3):
    tp = phase_triplet(mv_series[i], mv_series[j], mv_series[k])
    print(f'Triplet ({i+1},{j+1},{k+1}): {tp:.1f}°')
# All triplets are non-zero because moisture is monotonically decreasing
# (mv1 ≠ mv3 for all combinations)
```

---

## `fresnel_phase`

```python
fresnel_phase(mv, freq_ghz=1.4, theta_inc_deg=45.0, sand=0.51, clay=0.13)
```

### What it does

Computes the phase of the Fresnel transmission coefficient at the
air–soil interface, relative to dry soil (mv=0). This is the phase
jump the wave experiences when it crosses the boundary.

### The physics

When a wave crosses from air into soil, the boundary conditions require
the tangential electric and magnetic fields to be continuous. This gives
the Fresnel transmission coefficients:

**TE polarisation (HH — electric field horizontal):**
```
τ_TE = 2·kz_air / (kz_air + kz_soil')
```

**TM polarisation (VV — electric field in the incidence plane):**
```
τ_TM = 2·ε·kz_air / (ε·kz_air + kz_soil')
```

where:
- `kz_air = k₀·cos(θ_inc)` — vertical wavenumber in air
- `kz_soil'` — vertical wavenumber in soil (complex)
- `ε` — complex dielectric constant of soil

The **phase** of τ is the phase shift at the boundary. We return it
relative to dry soil (mv=0) so the plot shows only the moisture-induced
variation.

### Why this phase is small

The impedance contrast between air (ε=1) and soil (ε≈10) is large.
This means most of the wave is transmitted with almost no phase shift —
the boundary is nearly "transparent" in terms of phase.

The Fresnel phase variation with moisture is < 5°, compared to the
propagation phase (tens of degrees). This is why HH and VV
interferograms look almost identical for bare soil.

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv` | array_like | — | Volumetric soil moisture. |
| `freq_ghz` | float | `1.4` | Radar frequency. |
| `theta_inc_deg` | float | `45.0` | Incidence angle. |
| `sand` | float | `0.51` | Sand fraction. |
| `clay` | float | `0.13` | Clay fraction. |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `phase_te` | ndarray | TE (HH) transmission phase variation in degrees, relative to mv=0. |
| `phase_tm` | ndarray | TM (VV) transmission phase variation in degrees, relative to mv=0. |

### Example — reproduce De Zan (2014) Fig. 5

```python
import numpy as np
import matplotlib.pyplot as plt
from navasar.coherence import fresnel_phase

mv = np.linspace(0, 0.5, 200)
pte, ptm = fresnel_phase(mv, freq_ghz=1.4, theta_inc_deg=45.0)

plt.plot(mv, pte, 'b-', label='TE (HH)')
plt.plot(mv, ptm, 'r--', label='TM (VV)')
plt.xlabel('Volumetric moisture mv')
plt.ylabel('Transmission phase variation [deg]')
plt.title('Fresnel phase — De Zan (2014) Fig. 5')
plt.legend()
print(f'Max TE phase: {np.max(np.abs(pte)):.2f}°  ← negligible vs propagation phase')
```

---

## References

- De Zan, F., Parizzi, A., Prats-Iraola, P., & López-Dekker, P. (2014).
  *A SAR interferometric model for soil moisture.*
  IEEE Transactions on Geoscience and Remote Sensing, 52(1), 418–425.
  https://doi.org/10.1109/TGRS.2013.2241069

- Orfanidis, S. J. (2016). *Electromagnetic Waves and Antennas*, Chapter 7.
  http://www.ece.rutgers.edu/~orfanidi/ewa/
