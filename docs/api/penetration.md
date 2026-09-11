# `navasar.penetration` — Penetration Depth Models

Functions for computing one-way and two-way electromagnetic wave penetration depths in lossy dielectric soils.

---

## Functions

### `penetration_depth_dezan(mv, freq_ghz=1.4, inc_angle_deg=45.0, sand=0.51, clay=0.13)`
Computes penetration depth from vertical complex wavenumber $k'_z$ (De Zan et al. 2014, Eq. 5):
$$\delta = \frac{-1}{2\text{Im}(k'_z)}$$

### `penetration_depth_behari(mv, freq_ghz=1.4, sand=0.51, clay=0.13)`
Computes penetration depth from loss angle $\delta_L = \frac{1}{2}\arctan(\varepsilon'' / \varepsilon')$ (Behari 2005, Eq. 1.23):
$$P_d = \frac{\lambda}{4\pi \sqrt{|\varepsilon|} |\sin \delta_L|}$$

### `skin_depth_vs_moisture(mv_array, freq_ghz=1.4, method='dezan', **kwargs)`
Evaluates skin depth across an array of volumetric moisture contents.

### `penetration_depth_vs_frequency(mv, freq_array, **kwargs)`
Evaluates penetration depth across microwave frequency sweeps (e.g. 0.5 to 18 GHz).

### `two_way_penetration(mv, freq_ghz=1.4, **kwargs)`
Computes two-way round-trip radar penetration depth ($\delta / 2$).
