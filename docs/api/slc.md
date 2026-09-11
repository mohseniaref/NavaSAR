# `navasar.slc`

Single-look complex (SLC) SAR simulation.

**Source:** `navasar/slc.py`

---

## `simulate_slc`

```python
simulate_slc(mv_series, freq_ghz=1.4, sand=0.51, clay=0.13,
             n_pixels=400, n_layers=200, max_depth=0.30,
             sigma=None, rng=None)
```

Simulate a stack of SLC SAR acquisitions from a soil moisture time-series.

Zheng & Fattahi (2026), Eq. (3) and Table 1.

$$s_t = \sum_{n=0}^{N} \theta_n\,e^{-2j\,k^{(t)}_{z,n}\,z_n},
\quad \theta_n \sim \mathcal{CN}(0, \sigma_n),
\quad z_n = n\,D/N$$

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv_series` | array_like, shape `(T,)` or `(T, N)` | — | Volumetric soil moisture per acquisition. If 1-D, uniform moisture with depth is assumed. If 2-D, each row is the depth profile at that acquisition. |
| `freq_ghz` | float | `1.4` | Radar frequency in GHz |
| `sand` | float | `0.51` | Sand fraction |
| `clay` | float | `0.13` | Clay fraction |
| `n_pixels` | int | `400` | Number of pixels in the multi-look window ($M$) |
| `n_layers` | int | `200` | Number of depth layers ($N$) |
| `max_depth` | float | `0.30` | Maximum penetration depth in metres ($D$) |
| `sigma` | array_like or None | `None` | Radar cross-section per layer, shape `(N,)`. Uniform (all ones) if `None`. |
| `rng` | `np.random.Generator` or None | `None` | Random number generator. Uses `default_rng(42)` if `None`. |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `slc_stack` | complex ndarray, shape `(T, n_pixels)` | Simulated SLC stack. Each row is one acquisition. |

### Notes

- The scatterer realisations `theta` are drawn once and shared across all
  acquisitions, so temporal coherence is preserved.
- Increasing `n_pixels` reduces speckle noise in multi-looked products.
- Increasing `n_layers` improves the approximation of the continuous integral
  (Eq. 9 of De Zan 2014).
- For `n_pixels → ∞`, the multi-looked interferogram converges to the
  analytical expression of De Zan 2014, Eq. (12).

### Example

```python
import numpy as np
from navasar.slc import simulate_slc

mv_series = np.array([0.10, 0.20, 0.15, 0.25, 0.12])
slc_stack = simulate_slc(mv_series, freq_ghz=1.4, n_pixels=400, n_layers=200)
# slc_stack.shape == (5, 400)
```

### Varying moisture with depth

```python
T, N = 10, 200
# Desert: moisture increases linearly with depth
mv_depth = np.zeros((T, N))
for t in range(T):
    mv_depth[t] = np.linspace(0.05, 0.10, N)

slc_stack = simulate_slc(mv_depth, freq_ghz=1.4)
```

---

## `multilook_interferogram`

```python
multilook_interferogram(slc1, slc2)
```

Multi-looked interferogram from two SLC pixel arrays.

Zheng & Fattahi (2026), Eq. (4):

$$z_{12} = \frac{1}{M}\sum_{m=1}^{M} s_{1,m}\,s^*_{2,m}$$

### Parameters

| Name | Type | Description |
|------|------|-------------|
| `slc1` | complex ndarray, shape `(n_pixels,)` | SLC pixels at time $t_1$ |
| `slc2` | complex ndarray, shape `(n_pixels,)` | SLC pixels at time $t_2$ |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `z12` | complex scalar | Multi-looked interferogram |

### Notes

The normalised coherence magnitude is:

```python
coh = abs(z12) / np.sqrt(abs(multilook_interferogram(slc1, slc1))
                         * abs(multilook_interferogram(slc2, slc2)))
```

### Example

```python
from navasar.slc import simulate_slc, multilook_interferogram
import numpy as np

slc = simulate_slc(np.array([0.20, 0.25]), n_pixels=400)
z12 = multilook_interferogram(slc[0], slc[1])
z11 = multilook_interferogram(slc[0], slc[0])
z22 = multilook_interferogram(slc[1], slc[1])
coh = abs(z12) / np.sqrt(abs(z11) * abs(z22))
phase_deg = np.rad2deg(np.angle(z12))
```

---

## References

- Zheng, Y., & Fattahi, H. (2026).
  *Modeling, prediction, and retrieval of surface soil moisture from InSAR closure phase.*
  RSE, 333, 115104. https://doi.org/10.1016/j.rse.2025.115104

- De Zan, F., Parizzi, A., Prats-Iraola, P., & López-Dekker, P. (2014).
  *A SAR interferometric model for soil moisture.*
  IEEE TGRS, 52(1), 418–425. https://doi.org/10.1109/TGRS.2013.2241069
