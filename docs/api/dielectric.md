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
- Dry sandy loam: ε ≈ 2.9 + 0.4j
- Wet sandy loam (mv=0.3): ε ≈ 15 + 3j

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

and S = sand fraction, C = clay fraction (both in range 0–1).

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

**Real part ε' at C-band (5.0 GHz):**

| Coefficient | Constant term | × Sand (S) | × Clay (C) |
|-------------|--------------|------------|------------|
| A₀ | 2.927 | −0.012 | −0.001 |
| A₁ | 5.505 | +0.371 | +0.062 |
| A₂ | 114.826 | −0.389 | −0.547 |

**Imaginary part ε'' at C-band (5.0 GHz):**

| Coefficient | Constant term | × Sand (S) | × Clay (C) |
|-------------|--------------|------------|------------|
| B₀ | 0.004 | +0.001 | +0.002 |
| B₁ | 0.951 | +0.005 | −0.010 |
| B₂ | 15.354 | +0.143 | −0.261 |

For frequencies other than 1.4 or 5.0 GHz, the code linearly interpolates
between the two tables.

### Worked example by hand

Sandy loam: S = 0.51, C = 0.13, mv = 0.25, freq = 1.4 GHz.

Step 1 — compute A₀, A₁, A₂:
```
A₀ = 2.862 + (−0.012)×0.51 + 0.001×0.13 = 2.862 − 0.006 + 0.000 = 2.856
A₁ = 3.803 + 0.462×0.51 + (−0.341)×0.13 = 3.803 + 0.236 − 0.044 = 3.995
A₂ = 119.006 + (−0.500)×0.51 + 0.633×0.13 = 119.006 − 0.255 + 0.082 = 118.833
```

Step 2 — evaluate polynomial at mv = 0.25:
```
ε' = 2.856 + 3.995×0.25 + 118.833×0.0625
   = 2.856 + 0.999 + 7.427
   = 11.28
```

Step 3 — same for imaginary part:
```
B₀ = 0.356 − 0.003×0.51 − 0.008×0.13 = 0.354
B₁ = 5.507 + 0.044×0.51 − 0.002×0.13 = 5.529
B₂ = 17.753 − 0.313×0.51 + 0.206×0.13 = 17.620
ε'' = 0.354 + 5.529×0.25 + 17.620×0.0625 = 0.354 + 1.382 + 1.101 = 2.837
```

Result: ε ≈ 11.28 + 2.84j  ✓ (matches the code output)

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv` | array_like | — | Volumetric soil moisture. Range: 0 (bone dry) to 0.5 (saturated). Units: m³ water / m³ soil. |
| `freq_ghz` | float | `1.4` | Radar frequency in GHz. L-band = 1.4, C-band = 5.0. Other values are linearly interpolated. |
| `sand` | float | `0.51` | Sand fraction as a decimal (not percent). Default is the AGRISAR 2006 field 222 value. |
| `clay` | float | `0.13` | Clay fraction as a decimal. Note: sand + silt + clay = 1. |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `eps` | complex ndarray | ε = ε' + jε''. Real part controls phase velocity; imaginary part controls absorption. Both are positive. |

### Physical interpretation of the output

```python
eps = hallikainen_dielectric(0.25, freq_ghz=1.4)
# eps ≈ 11.28 + 2.84j

# Phase velocity in soil (fraction of speed of light):
v_phase = 1.0 / np.sqrt(eps.real)   # ≈ 0.30 → 30% of c

# Wavelength in soil (shorter than in air):
lam_soil = 0.3 / (1.4 * np.sqrt(eps.real))   # ≈ 6.4 cm vs 21.4 cm in air

# Loss tangent (ratio of absorption to storage):
tan_delta = eps.imag / eps.real   # ≈ 0.25 → significant loss
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
