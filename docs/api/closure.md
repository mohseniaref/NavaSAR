# `navasar.closure`

Closure phase time-series analysis.

**Source:** `navasar/closure.py`

---

## `closure_phase_timeseries`

```python
closure_phase_timeseries(mv_series, freq_ghz=1.4, sand=0.51, clay=0.13,
                          n_pixels=400, n_layers=200, max_depth=0.30,
                          rng=None)
```

Compute the closure phase time-series from a soil moisture time-series.

Zheng & Fattahi (2026): the closure phase is the difference between a
bandwidth-1 (nearest-neighbour) and a full-network InSAR time-series:

$$\text{CP}(t) = \phi^{\text{BW-1}}(t) - \phi^{\text{full-net}}(t)$$

### Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `mv_series` | array_like, shape `(T,)` or `(T, N)` | — | Volumetric soil moisture per acquisition. 1-D = uniform with depth; 2-D = depth profile per acquisition. |
| `freq_ghz` | float | `1.4` | Radar frequency in GHz |
| `sand` | float | `0.51` | Sand fraction |
| `clay` | float | `0.13` | Clay fraction |
| `n_pixels` | int | `400` | Multi-look window size ($M$) |
| `n_layers` | int | `200` | Number of depth layers ($N$) |
| `max_depth` | float | `0.30` | Maximum penetration depth in metres ($D$) |
| `rng` | `np.random.Generator` or None | `None` | Random number generator |

### Returns

| Name | Type | Description |
|------|------|-------------|
| `closure_ts` | ndarray, shape `(T,)` | Closure phase time-series in **degrees** |
| `bw1_ts` | ndarray, shape `(T,)` | BW-1 cumulative phase time-series in degrees |
| `fullnet_ts` | ndarray, shape `(T,)` | Full-network cumulative phase time-series in degrees |

### Notes

**BW-1 time-series** accumulates nearest-neighbour phases:

$$\phi^{\text{BW-1}}(t) = \sum_{i=0}^{t-1} \angle(z_{i,i+1})$$

**Full-network time-series** solves the least-squares system:

$$\mathbf{A}\,\boldsymbol{\phi} = \boldsymbol{\Delta\phi}$$

where $\mathbf{A}$ is the design matrix of all interferometric pairs and
$\boldsymbol{\Delta\phi}$ is the vector of multi-looked phases.

**Interpretation of closure phase signatures:**

| Soil moisture anomaly | Closure phase step |
|-----------------------|--------------------|
| Positive (rain event) | Positive step-change |
| Negative (drying) | Negative step-change |
| Duration of anomaly | Duration of step |
| Larger anomaly | Larger step amplitude |

### Example

```python
import numpy as np
from navasar.closure import closure_phase_timeseries

# Desert scenario: brief rain event
T = 20
mv_ts = np.full(T, 0.05)
mv_ts[6:10] = 0.18   # rain event

cp_ts, bw1_ts, fn_ts = closure_phase_timeseries(
    mv_ts, freq_ghz=1.4, n_pixels=400, n_layers=200
)

import matplotlib.pyplot as plt
plt.plot(cp_ts, label='Closure phase')
plt.axhline(0, color='k', ls='--')
plt.xlabel('Acquisition'); plt.ylabel('deg')
plt.legend(); plt.show()
```

### Varying moisture with depth

```python
T, N = 20, 200
mv_depth = np.zeros((T, N))
for t in range(T):
    base = 0.05 if t < 6 or t >= 10 else 0.18
    # Desert: moisture increases linearly with depth
    mv_depth[t] = np.linspace(base, base * 1.5, N)

cp_ts, _, _ = closure_phase_timeseries(mv_depth, freq_ghz=1.4)
```

---

## Internal functions

### `_bw1_timeseries(ifg_matrix)`

Cumulative phase from nearest-neighbour pairs only.

```
ifg_matrix : complex ndarray, shape (T, T)
    ifg_matrix[i, j] = multilook_interferogram(slc[i], slc[j])
Returns: ndarray, shape (T,)  — cumulative phase in radians
```

### `_fullnet_timeseries(ifg_matrix)`

Full-network time-series via least-squares inversion over all pairs.

```
ifg_matrix : complex ndarray, shape (T, T)
Returns: ndarray, shape (T,)  — cumulative phase in radians
```

---

## References

- Zheng, Y., & Fattahi, H. (2026).
  *Modeling, prediction, and retrieval of surface soil moisture from InSAR closure phase.*
  RSE, 333, 115104. https://doi.org/10.1016/j.rse.2025.115104

- Zheng, Y., et al. (2022).
  *Closure phase in InSAR processing: PSDS vs. SBAS.*
  IEEE TGRS. https://doi.org/10.1109/TGRS.2022.3167648
