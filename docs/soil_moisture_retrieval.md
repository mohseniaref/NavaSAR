# Soil-moisture retrieval: methods, physics, attribution and limits

**Executable tutorial:** [Notebook 15](../notebooks/15_soil_moisture_retrieval_comparison.ipynb).
Run all cells using the `sp` Python kernel. It saves timestamped PNG/PDF figures,
synthetic HDF5 arrays and provenance under `examples/15_soil_moisture_retrieval/`.
The optional real-data section reads existing point results only and does not
launch processing or change the source HDF5. Missing real data are skipped.

## Two different meanings of “model”

**Hallikainen, Dobson and Wang–Schmugge are dielectric models.** They predict
complex permittivity from moisture and ancillary soil properties. They are not,
by themselves, algorithms converting observed closure into moisture.

**Event SMI, gradient SMI and robust event SMI are retrieval algorithms.** The
first two are described in Zheng & Fattahi (2026); the third is our use of a
standard robust loss with the event model. Their normalized real-data results
do not presently use any dielectric model or per-pixel texture calibration.

| Method | Assumptions/input | Advantage | Limitation and attribution |
|---|---|---|---|
| Event SMI | Saturating exponential closure events | Smooth, interpretable event timing and decay | Published Eqs. 10–11; non-exponential histories can be poorly represented |
| Normalized gradient | Closure derivative after filtering | Fast; no explicit event-count fit | Published Section 6.1; sensitive to noise, sampling and smoothing |
| Robust event SMI | Same event model, soft-L1 loss | Can reduce influence of isolated outliers | Standard robust optimization; can also downweight genuine sharp events; no novelty claim |
| Dielectric closure inversion | Complex covariance, dielectric/scattering model, known moisture reference, soil parameters | Uses physical material response; works on wrapped triplets | Conditional and potentially nonunique; slower and sensitive to physical-model error; closure inversion predates this code |

## Is Event SMI exactly the paper's method?

The equations are the paper's main retrieval equations, implemented independently:

$$c(t)=c_0+\sum_i\beta_i(1-e^{-(t-t_i)/\omega_i})H(t-t_i),$$
$$x(t)=\sum_i\frac{\beta_i}{A\lambda\omega_i}
e^{-(t-t_i)/(\lambda\omega_i)}H(t-t_i).$$

With A = lambda = 1, following the paper, the code returns

$$\mathrm{SMI}(t)=\frac{x(t)-\min x}{\max x-\min x}.$$

This is **paper-equation-based**, not the authors' source code or an exact
reproduction of their entire workflow. Gradient-guided candidate search, BIC,
optimization bounds and robust loss are our numerical choices. Our pointwise
BW1-guided branch diagnostic also differs from the spatial phase-unwrapping and
references used in the paper's field workflow.

The paper's normalized-gradient proxy is

$$\mathrm{SMI}_g=\operatorname{normalize}_{[0,1]}
\left(\frac{d\,\operatorname{medianfilter}_3(c)}{dt}\right).$$

Actual elapsed days are used. A constant or missing signal is flagged, not
given artificial moisture. Normalization is not clipping, and a positive index
is not evidence of absolute moisture accuracy.

## NavaSAR APIs

```python
from navasar import closure_phase_from_ifgrams, retrieve_smi

# I: complete Hermitian complex multilooked pair-product matrix.
# phase sign: I_ij = mean(z_i * conjugate(z_j)).
closure_deg, bw1_deg, full_deg = closure_phase_from_ifgrams(I, branch="bw1")
standard = retrieve_smi(time_days, closure_deg)
robust = retrieve_smi(time_days, closure_deg, loss="soft_l1")
event_index = standard["event_smi"]
gradient_index = standard["gradient_smi"]
robust_index = robust["event_smi"]
```

`branch="bw1"` lifts to the nearest BW1-predicted branch; it is not spatial
unwrapping. Save and compare `branch="principal"` as a branch diagnostic.
Do not treat the difference as resolved simply because the arrays are finite.

`navasar.zheng_retrieval` also exposes `event_response` and `normalize_smi`.
SciPy is required for the retrieval modules. Core dielectric/SLC modules remain
importable without SciPy.

## Using dielectric models inside an inversion

The added `navasar.dielectric_inversion` module implements the following
**experimental alternative**, not an additional Zheng estimator:

$$\mathbf m \longrightarrow \epsilon(\mathbf m;\mathrm{texture},f,\ldots)
\longrightarrow k_i\longrightarrow C_{ij}\longrightarrow c_{ijk}.$$

It uses uniform semi-infinite volume scattering:

$$k_i=k_0\sqrt{\epsilon_i-\sin^2\theta},\quad\Im k_i<0,
\qquad C_{ij}=\frac{1}{2\mathrm{i}(k_i-k_j^*)}.$$

Normalize C to unit diagonal before computing

$$c_{ijk}=\arg(C_{ij}C_{jk}C_{ki}).$$

The code retains the local NavaSAR sign/root convention. This is not the
finite-depth multilayer simulator of Zheng–Fattahi. The assumed scattering
profile is just as important as the dielectric choice.

It minimizes weighted **circular** residuals rather than directly subtracting
unwrapped angles:

$$\min_{m_1,\ldots,m_{N-1}}\sum_q w_q
\left[\frac{\arg e^{\mathrm{i}(c_q^{\rm model}-c_q^{\rm obs})}}{0.1}\right]^2
+\eta^2\sum_i\left[\frac{m_i-m_{i-1}}
{0.1\sqrt{\Delta t_i/\operatorname{median}\Delta t}}\right]^2.$$

Acquisition 0 is fixed to the explicitly supplied `reference_moisture`. It must
come from independent evidence; do not silently assume 0.10 for real data.
The 0.1-radian residual scale and 0.1 moisture-increment scale are stated research
choices, not measured observation errors. Triplets are correlated. Multiple
starts, bounds, Jacobian rank and solution spread expose some ambiguities but
do not guarantee a unique solution. There is no held-out validation in this API.

```python
from navasar.dielectric_inversion import invert_closure_moisture

# Toy configuration at a supported frequency; NOT a Sentinel-1 calibration.
# reference_mv must be supplied from independent information.
result = invert_closure_moisture(
    I, time_days, reference_moisture=reference_mv,
    model="wang-schmugge",   # also "hallikainen" or "dobson"
    frequency_ghz=5.0, sand=0.51, clay=0.13,
    incidence_deg=40.0, bulk_density=1.3, temperature_c=20.0,
    bounds=(0.005, 0.45), smoothness=0.1, starts=3,
    allow_experimental=True,
)
conditional_moisture = result["moisture"]
print(result["model_note"], result["multistart_max_spread"])
```

Outputs have moisture-fraction coordinates **conditional on** the reference,
forward model, priors and assumed inputs. The returned flag
`absolute_moisture_validated` is deliberately false. Optimizer convergence is not
field validation. No existing real-result file is modified by importing/calling
these APIs. A later 2026-10-15 six-point direct-inversion experiment is described
in the linked source audit; it is separate from the saved real Event SMI.

### Why explicit experimental opt-in?

Source audit updated 2026-10-15 (see [detailed review](dielectric_model_review_20261015.md)):

- **Hallikainen:** corrected against MIMICS Table E.1: polynomial texture is
  percent by weight (API still accepts fractions). The old 5 GHz row was
  mislabeled 4 GHz; its loss coefficients were also wrong. Published rows
  1.4, 4, 6, 8, 10, 12, 14, 16, 18 GHz are now provided, with adjacent linear
  interpolation and no extrapolation. Sentinel-1 is now supported.
- **Dobson:** `behari.dobson_semiemp` remains the legacy texture-independent
  approximation. The inverse wrapper now uses `dobson_reviewed.dobson_ulaby`,
  the explicitly identified texture-aware Ulaby Module 4.7 textbook variant.
- **Wang–Schmugge:** texture is used, but density, temperature, conductivity and
  other empirical choices are not site measurements by default.

The 2026-10-15 update changes Hallikainen numerical results. Previously saved
synthetic figures/notebooks have not been overwritten and must be regenerated
before use with the reviewed coefficients. Existing real event indices do not
use these dielectric coefficients and are unchanged. A defensible physical
comparison requires verified formulations at the actual radar frequency,
consistent units/signs, a justified scattering profile, and sensitivity analysis
for soil/vegetation effects and the moisture reference.

## What is our contribution, and what is not?

### Displaying the real-point results

[Notebook 16](../notebooks/16_real_moisture_rainfall_maps.ipynb) separates direct
dielectric-model inversions from Event SMI, normalized gradient and robust Event
SMI. It displays rainfall overlays and geographic snapshots with common viridis
scales. The plotting entry point is
`MiaplPy/plot_moisture_methods_rain_maps.py`; it operates on saved H5 files and
does not perform a new inversion. Rainfall is comparison only, not an input or
prior. Six-point snapshots must not be interpreted as full-raster maps. No
spatial interpolation is performed. Relative indices are independently
normalized per point, so equal colors do not imply equal absolute moisture.

### Attribution

1. Event inversion is due to Zheng–Fattahi; the gradient proxy is also theirs.
2. Soft-L1 robust fitting is established numerical methodology. Its use here
   is an implementation extension, not evidence of a new scientific method.
3. Direct closure-based moisture inversion has prior work, including De Zan
   and Gomba (2018). Changing dielectric models does not establish novelty.
4. Our current contribution is software integration, explicit provenance,
   controlled comparisons, geometry-aware sampling and diagnostics.

We cannot establish that nobody has used this particular combination before.
A publication-level improvement would require a focused literature review and
repeated synthetic and independent field evaluation demonstrating an advantage.
In the completed single known-event synthetic test, the gradient proxy had
lower normalized error than either event fit; robust fitting was not generally
demonstrated superior. The local dielectric forward tests gave poor recovery.

## Tests and research provenance

`tests/test_soil_retrieval.py` checks event recovery, zero-information behavior,
covariance Hermitian/PSD properties, experimental gating and conditional fits
for each local dielectric model. The eight-date noise-free inverse tests use
12 starting points: two starts initially gave 4.6–6.5 degree closure errors,
whereas 12 starts reduced those errors below 0.004 degrees. This demonstrates
local-minimum sensitivity, not a guarantee that 12 starts suffice on field data.
Numerical self-consistency is not physical
validation. The more extensive application tests and report are in
`../MiaplPy/test_zheng_paper_comparison.py` and
`../MiaplPy/docs/zheng_paper_comparison.md` relative to the repository root's parent.

Run from the NavaSAR repository with SciPy installed:

```bash
PYTHONPATH=. python -m unittest discover -s tests -p test_soil_retrieval.py -v
```

## Primary sources

- [Zheng & Fattahi (2026), Eqs. 10–11 and Section 6.1](https://doi.org/10.1016/j.rse.2025.115104).
- [De Zan & Gomba (2018), closure-phase inversion](https://elib.dlr.de/121676/).
- [De Zan & Gomba (2018), inversion ambiguities](https://elib.dlr.de/120480/1/dezan_eusar_2018_closure_phase_ambiguities.pdf).
- [SciPy least_squares: established robust losses](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html).
