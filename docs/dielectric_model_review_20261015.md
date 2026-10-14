# Soil-moisture model audit and comparison — 15 October 2026

This audit is included with NavaSAR. The paper and result files referenced below
belong to the sibling local `MiaplPy` project and are not distributed here.
On GitHub, use notebooks 15–16 for embedded examples and their limitations;
local result links require the original research workspace.

## What the local Zheng–Fattahi paper actually does

Source: `paper/1-s2.0-S0034425725005085-main.pdf`, Sections 2–4, 6.1–6.3,
Appendices A and C. Appendix A explicitly uses Hallikainen **1.4, 4, 6,
10 GHz** for L/S/C/X bands. Its C-band simulations use 6 GHz, acknowledging
the difference from Sentinel-1's 5.4 GHz. Saying that the published model
stops at 5 GHz was incorrect.

The paper retrieves an entire time series, not just selected dates. It fits
closure changes using a sum of event responses (Eq. 10), then constructs the
corresponding moisture precursor (Eq. 11) and normalizes it. For elapsed time t:

$$c(t)=c_0+\sum_k\beta_k[1-\exp(-(t-t_k)/\omega_k)]H(t-t_k),$$

$$x(t)=\sum_k\frac{\beta_k}{\omega_k}\exp(-(t-t_k)/\omega_k)H(t-t_k),
\qquad \mathrm{SMI}(t)=\frac{x(t)-\min x}{\max x-\min x}.$$

The second expression uses the paper's simplifying A=lambda=1 convention.
Unknown amplitude and time-warp parameters prevent interpreting this index
as absolute m³/m³. The gradient proxy is also discussed in the paper; it is
not a new contribution. A zero dynamic range produces an undefined index,
not a fabricated zero-moisture curve.

Our event-response equations follow the paper, but the optimizer, BIC event
selection and BW1-guided phase-branch lifting are implementation choices.
The existing six-point experiment is **not an exact end-to-end reproduction**
of the paper's spatially unwrapped, referenced and spatially averaged field
workflow. The large principal/BW1-branch differences remain a diagnostic issue.

## Three dielectric models are forward-model choices, not three SMI algorithms

They map moisture/texture/physical inputs to complex permittivity. To retrieve
moisture directly, they must be combined with a scattering model and an inverse
problem. Changing dielectric models is therefore meaningful for physical
closure inversion and simulations, but does not automatically change the
paper's empirical event-SMI fit to the same measured closure series.

### Hallikainen: verified corrections implemented

[Original article, DOI](https://doi.org/10.1109/TGRS.1985.289497).
Coefficient transcription checked visually against the author-group
[MIMICS report](https://web.eecs.umich.edu/~lep/mimics.pdf), printed pages
E-3/E-4, PDF pages 146/147, Eq. E.14 and Table E.1.

For each real/loss component:

$$\epsilon=(a_0+a_1S+a_2C)+(b_0+b_1S+b_2C)m_v
+(c_0+c_1S+c_2C)m_v^2.$$

S and C are sand/clay **percent by weight**, not fractions. Our public API
continues to accept fractions and now converts them internally. The previous
row labelled 5 GHz used 4 GHz real coefficients; its quadratic loss row also
did not match the source. Both are corrected. All nine published table rows
from 1.4 to 18 GHz are implemented. At 5.405 GHz we interpolate adjacent 4/6
GHz predictions; this interpolation is our numerical choice, not a new
measured coefficient row. No out-of-range extrapolation is allowed.

Hand-calculated regression checks at m=0.2, sand=0.51, clay=0.13:

| Frequency | Real permittivity | Loss component |
|---|---:|---:|
| 4 GHz | 10.86344 | 1.50904 |
| 6 GHz | 10.24868 | 1.94732 |

NavaSAR uses epsilon'=real plus i*loss; the report writes minus i*loss.
The existing propagation/root convention is retained consistently. Empirical
polynomials can predict nonphysical loss for some dry/extreme-texture inputs;
the inverse wrapper rejects non-attenuating responses, rather than flooring them.

### Dobson: texture-aware textbook variant added

[Original mixing-model article](https://doi.org/10.1109/TGRS.1985.289498).
The old `behari.dobson_semiemp` uses a fixed exponent and ignores texture.
It remains available for historical reproduction, but the inverse wrapper
now uses `dobson_reviewed.dobson_ulaby` from the primary textbook companion
[Module 4.7 implementation](https://ulabyweb.eecs.umich.edu/mrs/codes/Module4_7/plotLibrary.js).

This specific high-frequency variant uses fraction texture S,C:

$$\beta_1=1.27-0.519S-0.152C,\quad
\beta_2=2.06-0.928S-0.255C,\quad\alpha=0.65,$$

$$\sigma=-1.645+1.939\rho_b-2.256S+1.594C,$$

$$\epsilon'=\left[1+0.66\rho_b+m^{\beta_1}(\epsilon'_w)^\alpha-m\right]^{1/\alpha},
\qquad\epsilon''=m^{\beta_2}\epsilon''_w.$$

The water loss includes conductivity scaled by porosity/m, frequency and vacuum
permittivity. The dry limit is evaluated without division by zero. Temperature
and density must be supplied or clearly identified as assumptions. We implement
only 1.4–18 GHz and 0–40 °C here; this does not claim equivalence to every
Dobson/Peplinski variant. Negative modeled loss and moisture above porosity are
rejected. Legacy results and this variant should not be mixed without labels.

### Wang–Schmugge: useful but not parameter-free

The [original NASA report](https://ntrs.nasa.gov/api/citations/19790010347/downloads/19790010347.pdf)
explicitly models a bound/free-water transition linked to texture. Its wilting
point equation uses percent-weight sand/clay. Density controls porosity; water
permittivity depends on temperature/frequency. Knowing sand/clay alone does not
determine all these quantities.

Our current comparison retains the local empirical implementation, with assumed
1.3 g/cm³ and 20 °C. Its precise loss choices have **not** completed the same
source-by-source audit as Hallikainen. It must be labelled a local variant, not
an exact original-paper reproduction. No new site density/temperature data
were downloaded in this task.

## Alternative inversion experiment

Section 6.3 of Zheng–Fattahi discusses prior physical inversion and empirical
approaches. [De Zan–Gomba (2018)](https://doi.org/10.1016/j.rse.2018.08.034)
predates our physical closure inversion. [Wig et al. (2024)](https://doi.org/10.1109/TGRS.2024.3399069)
also predates this work. We have not implemented an exact Wig retrieval here.

The new comparison script fits wrapped triangle closures for all 79 dates at
six saved points, using three dielectric choices and a uniform semi-infinite
scattering covariance. It uses multistart, bounds and time-scaled smoothing.
An initial moisture of 0.10 m³/m³ is explicitly **assumed**, not measured or
estimated from ERA5. Results are conditional on that anchor and scattering
assumption. Closure residual, bound fraction, Jacobian rank, optimizer status
and multistart spread are preserved. This is a sensitivity experiment, not
evidence of novel or uniquely identified absolute soil moisture.

Completed results: [three-model report](../results/dielectric_point_comparison_20261015_reviewed/report.md)
and [figure](../results/dielectric_point_comparison_20261015_reviewed/fig01_conditional_three_model_timeseries.png).
Of 18 attempted point/model fits, six converged, eleven reached the iteration
limit, and Dobson at west_3 was rejected because a trial moisture yielded
negative empirical dielectric loss. Multistart spreads were approximately
0.345–0.445 m³/m³. These are different optimizer solutions, not uncertainty
confidence intervals. No physical winner is established and full-raster
inversion with these settings is not justified by this test.

## ERA5-Land comparison completed

See [results and correlation table](../results/soil_moisture_era5land_20261015/report.md).
All 79 dates and six points were compared with surface 0–7 cm ERA5-Land;
7–28 cm and daily means are additional checks. Hourly sampling uses 09:44:53 UTC
from the user-supplied original SLC metadata, assumed constant across dates.
No InSAR-to-ERA5 calibration is performed.

The largest raw Event-SMI correlation is east_4, rho=0.460, but its
calendar-month-demeaned correlation is 0.023. West_1 Event-SMI has raw 0.228
and month-demeaned 0.262. Results are mixed, not validation of a winning model.
Four ERA5-Land grid cells represent six points: east_5/east_6 share one and
west_2/west_3 another. Pixel-scale differences cannot be resolved by these
coarse reanalysis cells. Normalized SMI and volumetric moisture have different
units, so no misleading absolute RMSE is reported.

## Reproducibility and previous outputs

Code: `compare_era5land_moisture.py`, `compare_dielectric_point_inversions.py`.
ERA5 raw JSON, request URLs, actual grid coordinates, H5 time series, PNG/PDF
figures and Markdown reports are retained. Direct-inversion results record
source hashes. Existing results are not overwritten by default.

Historical Hallikainen-based simulation/notebook figures predate the coefficient
correction and require regeneration. The saved real event indices are derived
from observed covariances and are not modified by this correction. This task
does not launch a full-raster inversion or label unfinished numerical audits
as publication-ready validation.
