# `navasar.inversion` — Soil Moisture Retrieval & Inversion

Routines for inverting volumetric soil moisture from interferometric coherence and closure phase observables.

**Current retrieval guide:** [methods, actual APIs, attribution and limitations](../soil_moisture_retrieval.md).
Use `navasar.zheng_retrieval.retrieve_smi` for paper-equation event SMI and the
normalized-gradient proxy; `loss="soft_l1"` selects the robust extension.
Use `navasar.dielectric_inversion.invert_closure_moisture` for the separate
experimental material-model inverse problem. The legacy routines below are not
the current event/normalized-gradient SMI API. Signatures list key arguments;
consult Python help for additional simulation parameters.

---

## Coherence Inversion (De Zan 2014, Sec. III-D)

### `calibration_coherence_curve(mv_ref, freq_ghz=1.4, theta_inc_deg=45.0, sand=0.51, clay=0.13, n_points=100)`
Precomputes lookup tables of $|\gamma|(m_{v2})$ and $\angle\gamma(m_{v2})$ for a given reference master moisture.

### `invert_moisture_from_coherence(gamma_obs, mv_ref, freq_ghz=1.4, theta_inc_deg=45.0, sand=0.51, clay=0.13, observable='magnitude')`
Inverts slave soil moisture $m_{v2}$ via 1D Brent root-finding against the forward model.

---

## Historical transfer experiments — not the paper's main event inversion

The rectangular calibration below may have zero closure and now raises rather
than manufacturing a transfer. It must not be used as absolute calibration.

### `calibrate_closure_transfer(mv_base, freq_ghz=1.4, ..., anomaly_amps=None, tau_acq=4, rng=None)`
Historical rectangular-event experiment attempting to fit a transfer coefficient
$T_{\text{coeff}}$ [deg / ($m^3/m^3$)] relating closure phase steps to anomalies:
$$\beta \approx T_{\text{coeff}} \times \alpha$$

### `insar_smi_timeseries(cp_ts_deg, T_coeff, smooth_window=3)`
Legacy acquisition-index gradient divided by a supplied transfer; this is not
the normalized-gradient proxy API and does not itself provide calibration:
$$\text{SMI}(t) = \frac{\nabla_t \phi_{\text{closure}}(t)}{T_{\text{coeff}}}$$

### `invert_moisture_from_triplet(triplet_step_deg, mv_base, freq_ghz=1.4, ...)`
Legacy step/transfer division using `calibrate_closure_transfer`; not a
three-date nonlinear dielectric inversion. Zero transfer now raises an error.

### `joint_inversion(gamma_mag_obs, triplet_obs_deg, mv_ref, freq_ghz=1.4, theta_inc_deg=45.0, sand=0.51, clay=0.13, weights=(1.0, 1.0))`
Retrieves slave moisture by minimizing a combined weighted magnitude and phase loss function.
