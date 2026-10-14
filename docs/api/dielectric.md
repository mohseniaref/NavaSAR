# `navasar.dielectric` — Soil Dielectric Properties

**Source:** `navasar/dielectric.py`

This module answers one question: **given a soil moisture value, what is
the dielectric constant of the soil, and how does the radar wave travel
through it?**

---

## Background: why does moisture change the dielectric constant?

Imagine filling a bucket with dry sand (dielectric constant ε ≈ 4).
Now slowly add water (ε ≈ 80). The mixture's ε rises because you are
replacing air (ε = 1) with a material that is 80× more polarisable.

The dielectric constant ε is a **complex number**:

```
ε = ε' + j·ε''
     ↑       ↑
  real part  imaginary part
  (slows     (absorbs
   the wave)  the wave)
```

Both parts increase as the soil gets wetter. At L-band (1.4 GHz):
- Dry soil (51% sand, 13% clay): ε ≈ 2.26 + 0.10j
- Wet soil with the same texture (mv=0.3): ε ≈ 18.30 + 2.82j

---

## `hallikainen_dielectric`

```python
hallikainen_dielectric(mv, freq_ghz=1.4, sand=0.51, clay=0.13)
```

### What it does

Returns the complex dielectric constant of wet soil using the empirical
polynomial model of Hallikainen et al. (1985). This is the standard
model used in De Zan (2014) and Zheng & Fattahi (2026).

### The formula

The model fits a **second-order polynomial in mv** separately for the
real and imaginary parts:

```
ε'(mv) = A₀ + A₁·mv + A₂·mv²
ε''(mv) = B₀ + B₁·mv + B₂·mv²
```

where each coefficient itself depends on soil texture:

```
A₀ = a₀₀ + a₀₁·S + a₀₂·C
A₁ = a₁₀ + a₁₁·S + a₁₂·C
A₂ = a₂₀ + a₂₁·S + a₂₂·C
```

and S = sand percentage, C = clay percentage (by weight, 0–100).
The API accepts fractions and multiplies by 100 internally.

### Full coefficient tables

**Real part ε' at L-band (1.4 GHz):**

| Coefficient | Constant term | × Sand (S) | × Clay (C) |
|-------------|--------------|------------|------------|
| A₀ (offset) | 2.862 | −0.012 | +0.001 |
| A₁ (× mv) | 3.803 | +0.462 | −0.341 |
| A₂ (× mv²) | 119.006 | −0.500 | +0.633 |

**Imaginary part ε'' at L-band (1.4 GHz):**

| Coefficient | Constant term | × Sand (S) | × Clay (C) |
|-------------|--------------|------------|------------|
| B₀ (offset) | 0.356 | −0.003 | −0.008 |
| B₁ (× mv) | 5.507 | +0.044 | −0.002 |
| B₂ (× mv²) | 17.753 | −0.313 | +0.206 |

**Real part ε' at 4.0 GHz (previously mislabeled 5.0):**

| Coefficient | Constant term | × Sand (S) | × Clay (C) |
|-------------|--------------|------------|------------|
| A₀ | 2.927 | −0.012 | −0.001 |
| A₁ | 5.505 | +0.371 | +0.062 |
| A₂ | 114.826 | −0.389 | −0.547 |

**Imaginary part ε'' at 4.0 GHz:**

| Coefficient | Constant term | × Sand (S) | × Clay (C) |
|-------------|--------------|------------|------------|
| B₀ | 0.004 | +0.001 | +0.002 |
| B₁ | 0.951 | +0.005 | −0.010 |
| B₂ | 16.759 | +0.192 | +0.290 |

The implementation includes 1.4, 4, 6, 8, 10, 12, 14, 16 and 18 GHz rows.
Other frequencies within this range interpolate adjacent rows (5.405 uses 4/6).
No extrapolation. Tables checked against [MIMICS E.1](https://web.eecs.umich.edu/~lep/mimics.pdf).
This correction changes historical simulations; previous figures need rerunning.

### Worked example by hand

API sand=0.51, clay=0.13 means polynomial S=51, C=13; mv=0.25, freq=1.4 GHz.

Step 1 — compute A₀, A₁, A₂:
```
A₀ = 2.862 − 0.012×51 + 0.001×13 = 2.263
A₁ = 3.803 + 0.462×51 − 0.341×13 = 22.932
A₂ = 119.006 − 0.500×51 + 0.633×13 = 101.735
```

Step 2 — evaluate polynomial at mv = 0.25:
```
ε' = 2.263 + 22.932×0.25 + 101.735×0.0625 = 14.35444
```

Step 3 — same for imaginary part:
```
B₀ = 0.356 − 0.003×51 − 0.008×13 = 0.099
B₁ = 5.507 + 0.044×51 − 0.002×13 = 7.725
B₂ = 17.753 − 0.313×51 + 0.206×13 = 4.468
ε'' = 0.099 + 7.725×0.25 + 4.468×0.0625 = 2.3095
```

Result: ε ≈ 14.35 + 2.31j.

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv` | array_like | — | Volumetric soil moisture. Range: 0 (bone dry) to 0.5 (saturated). Units: m³ water / m³ soil. |
| `freq_ghz` | float | `1.4` | Radar frequency in GHz, 1.4–18. Adjacent-table interpolation. |
| `sand` | float | `0.51` | Sand fraction as a decimal (not percent). Default is the AGRISAR 2006 field 222 value. |
| `clay` | float | `0.13` | Clay fraction as a decimal. Note: sand + silt + clay = 1. |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `eps` | complex ndarray | ε = ε' + jε''. Empirical fits can become nonphysical at extreme inputs; the inversion rejects non-attenuating results. |

### Physical interpretation of the output

```python
eps = hallikainen_dielectric(0.25, freq_ghz=1.4)
# eps ≈ 14.35 + 2.31j

# Phase velocity in soil (fraction of speed of light):
v_phase = 1.0 / np.sqrt(eps.real)   # ≈ 0.264 (low-loss approximation)

# Wavelength in soil (shorter than in air):
lam_soil = 0.3 / (1.4 * np.sqrt(eps.real))   # ≈ 5.7 cm vs 21.4 cm in air

# Loss tangent (ratio of absorption to storage):
tan_delta = eps.imag / eps.real   # ≈ 0.161
```

### Example

```python
import numpy as np
import matplotlib.pyplot as plt
from navasar.dielectric import hallikainen_dielectric

mv = np.linspace(0, 0.5, 200)
eps = hallikainen_dielectric(mv, freq_ghz=1.4, sand=0.51, clay=0.13)

plt.plot(mv, eps.real, label="ε' (real)")
plt.plot(mv, eps.imag, label="ε'' (imaginary)")
plt.xlabel('Volumetric moisture mv')
plt.ylabel('Dielectric constant')
plt.legend()
plt.title('Hallikainen 1985 — L-band, sandy loam')
# Reproduces De Zan (2014) Fig. 2
```

---

## `kz_soil`

```python
kz_soil(mv, freq_ghz=1.4, theta_inc_deg=45.0, sand=0.51, clay=0.13)
```

### What it does

Computes the **vertical complex wavenumber** $k'_z$ inside the soil.
This is the key quantity that links the dielectric constant to the
interferometric phase and coherence.

### Step-by-step derivation

**Step 1 — free-space wavenumber:**

The wavenumber in air is:
```
k₀ = 2π / λ = 2π × f / c
```
At 1.4 GHz: k₀ = 2π × 1.4×10⁹ / 3×10⁸ = 29.3 rad/m

**Step 2 — horizontal wavenumber (conserved across the boundary):**

Snell's law says the horizontal component of k does not change when
the wave crosses the air–soil interface:
```
kx = k₀ × sin(θ_inc)
```
At θ_inc = 45°: kx = 29.3 × sin(45°) = 20.7 rad/m

**Step 3 — vertical wavenumber in soil:**

Inside the soil, the wave equation requires:
```
kx² + kz'² = k₀² × ε
```
Solving for kz':
```
kz' = sqrt(k₀² × ε − kx²)
```
Because ε is complex, kz' is also complex.

**Step 4 — choose the physical root:**

The square root has two solutions. We need the one where the wave
**attenuates downward** (not grows). This means Im(kz') must be negative:
```
if Im(kz') > 0:  flip sign → kz' = −kz'
```
This ensures the wave amplitude e^(−j·kz'·z) decays as z increases.

### What the output tells you

```
kz' = Re(kz') + j·Im(kz')
       ↑              ↑
  phase per metre   attenuation per metre
  (controls         (Im < 0, so amplitude
   interferometric   decays with depth)
   phase)
```

**Penetration depth** (depth where amplitude falls to 1/e):
```
δ = −1 / (2 × Im(kz'))
```

**Two-way phase accumulated to depth z:**
```
φ(z) = −2 × Re(kz') × z   [radians]
```
The factor of 2 is because the wave travels down to depth z and back up.

### Numerical example

At mv=0.25, freq=1.4 GHz, θ=45°, sandy loam:
```
ε ≈ 11.28 + 2.84j
k₀ = 29.3 rad/m
kx = 20.7 rad/m
k₀²·ε = 29.3² × (11.28 + 2.84j) = 9672 + 2436j
kz'² = 9672 + 2436j − 20.7² = 9672 + 2436j − 428 = 9244 + 2436j
kz' ≈ −97.1 − 12.6j  rad/m

Penetration depth δ = −1/(2×(−12.6)) = 4.0 cm
Phase per cm of depth = 97.1 × 0.01 = 0.97 rad = 55.5°
```

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv` | array_like | — | Volumetric soil moisture |
| `freq_ghz` | float | `1.4` | Radar frequency in GHz |
| `theta_inc_deg` | float | `45.0` | Incidence angle in degrees from vertical. Typical SAR: 20°–50°. |
| `sand` | float | `0.51` | Sand fraction |
| `clay` | float | `0.13` | Clay fraction |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `kz` | complex ndarray | Vertical wavenumber [rad/m]. Re(kz) > 0 controls phase; Im(kz) < 0 controls attenuation. |

### Example

```python
import numpy as np
from navasar.dielectric import kz_soil

mv = np.linspace(0.01, 0.50, 100)
kz = kz_soil(mv, freq_ghz=1.4, theta_inc_deg=45.0)

# Penetration depth in cm
depth_cm = -100.0 / (2 * kz.imag)

# Phase accumulated over 5 cm depth (one-way), in degrees
phase_5cm_deg = np.rad2deg(-kz.real * 0.05)

print(f'At mv=0.10: depth={depth_cm[10]:.1f} cm, phase/5cm={phase_5cm_deg[10]:.1f}°')
print(f'At mv=0.30: depth={depth_cm[60]:.1f} cm, phase/5cm={phase_5cm_deg[60]:.1f}°')
```

---

## References

- Hallikainen, M. T., Ulaby, F. T., Dobson, M. C., El-Rayes, M. A., & Wu, L.-K. (1985).
  *Microwave dielectric behavior of wet soil – Part 1: Empirical models and experimental observations.*
  IEEE Transactions on Geoscience and Remote Sensing, GE-23(1), 25–34.
  https://doi.org/10.1109/TGRS.1985.289497

- De Zan, F., Parizzi, A., Prats-Iraola, P., & López-Dekker, P. (2014).
  *A SAR interferometric model for soil moisture.*
  IEEE Transactions on Geoscience and Remote Sensing, 52(1), 418–425.
  https://doi.org/10.1109/TGRS.2013.2241069

- Orfanidis, S. J. (2016). *Electromagnetic Waves and Antennas*, Chapter 7.
  http://www.ece.rutgers.edu/~orfanidi/ewa/
