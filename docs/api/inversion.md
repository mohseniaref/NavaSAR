# `navasar.inversion` — Soil Moisture Retrieval & Inversion

Routines for inverting volumetric soil moisture from interferometric coherence and closure phase observables.

---

## Coherence Inversion (De Zan 2014, Sec. III-D)

### `calibration_coherence_curve(mv_master, freq_ghz=1.4, n_points=200, sand=0.51, clay=0.13)`
Precomputes lookup tables of $|\gamma|(m_{v2})$ and $\angle\gamma(m_{v2})$ for a given reference master moisture.

### `invert_moisture_from_coherence(obs_value, mv_master, freq_ghz=1.4, target='magnitude', sand=0.51, clay=0.13)`
Inverts slave soil moisture $m_{v2}$ via 1D Brent root-finding against the forward model.

---

## Closure Phase InSAR-SMI (Zheng & Fattahi 2026, Sec. 4)

### `calibrate_closure_transfer(mv_base=0.05, anomaly_range=(0.02, 0.25), n_points=8, freq_ghz=1.4, **sim_kw)`
Calibrates the linear transfer coefficient $T_{\text{coeff}}$ [deg / ($m^3/m^3$)] relating closure phase steps to moisture anomalies:
$$\beta \approx T_{\text{coeff}} \times \alpha$$

### `insar_smi_timeseries(closure_phase_ts, transfer_coeff)`
Computes the InSAR Soil Moisture Index (SMI) time series from the temporal gradient of the closure phase:
$$\text{SMI}(t) = \frac{\nabla_t \phi_{\text{closure}}(t)}{T_{\text{coeff}}}$$

### `invert_moisture_from_triplet(triplet_deg, mv1, mv3, freq_ghz=1.4, sand=0.51, clay=0.13)`
Retrieves intermediate acquisition moisture $m_{v2}$ from a single phase triplet $\phi_{123}$.

### `joint_inversion(obs_coherence, obs_phase_deg, mv_master, freq_ghz=1.4, weight_mag=1.0, weight_phase=0.1, sand=0.51, clay=0.13)`
Retrieves slave moisture by minimizing a combined weighted magnitude and phase loss function.
