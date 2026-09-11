# Physical Theory — Comprehensive Primer & Reference

This document provides a thorough, step-by-step physical and mathematical derivation of all models implemented in `navasar`.
It is structured to guide undergraduates, graduate students, and researchers alike through the entire chain of physics:
from individual molecular dipoles, to heterogeneous soil mixture electrodynamics, to SAR volume scattering, and finally to interferometric closure phase retrieval algorithms.

---

## 0. The Big Picture

Synthetic Aperture Radar (SAR) systems emit coherent microwave pulses toward Earth. When a radar wave reaches bare ground, it doesn't reflect merely off an infinitesimally thin optical boundary. Instead, the electromagnetic radiation **penetrates several centimetres into the soil volume** before scattering back to the sensor.

### Key Mechanism
1. **Dry soil particles** (silicates, quartz, clay lattices) have a relatively low dielectric constant: $\varepsilon_r \approx 3 - 4$.
2. **Liquid water** consists of polar molecules that rotate under alternating electric fields, resulting in a huge dielectric permittivity: $\varepsilon_r \approx 80$ at room temperature and microwave frequencies.
3. Adding even a small fraction of water (e.g. 5% to 20% volumetric moisture $m_v$) dramatically raises the soil bulk permittivity.
4. Higher permittivity slows down wave velocity ($v = c / \sqrt{\varepsilon_r}$) and increases attenuation.
5. In SAR interferometry (comparing two radar passes), this changes both the **interferometric coherence magnitude** $|\gamma|$ (decorrelation) and the **interferometric phase** $\angle\gamma$ (propagation delay).
6. Because moisture-induced phase shifts behave non-linearly with respect to moisture content, three-pass interferograms yield a **non-zero closure phase** (phase triplet $\phi_{123} \neq 0$). This non-closure signature serves as an unambiguous tracer for soil moisture dynamics!

---

## 1. Dielectric Permittivity: Physics and Notation

Every material's response to an applied electromagnetic field $\mathbf{E}$ is governed by its complex relative permittivity:

$$\varepsilon = \varepsilon' + j \varepsilon''$$

*(Note: in engineering and physics literature, either $+j\omega t$ with $\varepsilon = \varepsilon' - j\varepsilon''$ or $-j\omega t$ is adopted. In `navasar` and De Zan et al., the time convention is $\exp(-j\omega t)$, meaning loss is positive imaginary or handled explicitly with negative imaginary spatial wavenumbers).*

- **Real part $\varepsilon'$ (Permittivity / Storage):** Measures polarization and the degree to which the medium slows down the wave.
- **Imaginary part $\varepsilon''$ (Dielectric Loss):** Measures energy dissipation (Ohmic absorption and molecular relaxation damping).
- **Loss Tangent:** $\tan \delta = \varepsilon'' / \varepsilon'$.

---

## 2. Hallikainen et al. (1985) Empirical Soil Model

Hallikainen et al. measured soil samples across 1.4 to 18 GHz and established empirical polynomials relating soil texture and moisture $m_v$ to permittivity:

$$\varepsilon'_r = (a_0 + a_1 S + a_2 C) + (b_0 + b_1 S + b_2 C) m_v + (c_0 + c_1 S + c_2 C) m_v^2$$
$$\varepsilon''_r = (a_0'' + a_1'' S + a_2'' C) + (b_0'' + b_1'' S + b_2'' C) m_v + (c_0'' + c_1'' S + c_2'' C) m_v^2$$

### Variable & Symbol Glossary
- $m_v$: Volumetric soil moisture [$m^3/m^3$], typical range $0.02 - 0.45$.
- $S$: Sand mass fraction [$0 - 1$] (grain size $0.05 - 2.0$ mm).
- $C$: Clay mass fraction [$0 - 1$] (grain size $< 0.002$ mm).
- $a_k, b_k, c_k$: Empirically calibrated polynomial coefficients (Table 1 of Hallikainen 1985).

### Why Texture Matters: Bound vs. Free Water
Clay platelets possess enormous specific surface areas ($~800\text{ m}^2/\text{g}$) with negative electrical surface charges. Water molecules within the first 2–3 monolayers are tightly electrostatically bound ("bound water"), preventing their free rotation. Thus, clay-rich soils exhibit lower permittivity than sandy soils at identical volumetric moisture contents!

---

## 3. Wave Propagation in Soil (De Zan 2014, Eq. 5)

When a plane radar wave in air ($k_0 = 2\pi / \lambda = \omega/c$) impinges upon a flat soil interface at incidence angle $\theta_{inc}$, boundary continuity dictates conservation of horizontal wavenumber:

$$k_x = k_0 \sin \theta_{inc}$$

From the dispersion relation $k_x^2 + k_z^2 = k_0^2 \varepsilon_r$, the vertical wavenumber inside the soil is:

$$k'_z(\varepsilon) = \sqrt{k_0^2 \varepsilon_r - k_x^2}$$

### Physical Root Selection & Attenuation
Since $\varepsilon_r$ has an imaginary loss component, $k'_z$ is complex:

$$k'_z = \text{Re}(k'_z) + j \text{Im}(k'_z)$$

To guarantee that downward-propagating waves $\exp(-j k'_z z)$ attenuate rather than amplify with depth $z$, **the physical branch requires $\text{Im}(k'_z) < 0$**.

---

## 4. Penetration Depth & Surface Emissivity (Behari 2005, Ch. 1)

### One-Way Penetration Depth (De Zan vs. Behari)
1. **De Zan Wavenumber Formulation:**
   The depth at which power drops to $1/e$ (amplitude drops to $e^{-1/2}$):
   $$\delta = \frac{-1}{2 \text{Im}(k'_z)}$$
2. **Behari Formulation (Eq. 1.23):**
   Expressed in terms of the dielectric loss angle $\delta_L = \frac{1}{2} \arctan(\varepsilon'' / \varepsilon')$:
   $$P_d = \frac{\lambda}{4\pi \sqrt{|\varepsilon|} |\sin \delta_L|}$$

### Surface Emissivity (Behari Eq. 1.27, 1.30)
By Kirchhoff's law of thermal radiation, emissivity $e$ and power reflectivity $r$ satisfy $e + r = 1$.
- **Smooth Surface:**
  $$e_s = 1 - \left| \frac{1 - \sqrt{\varepsilon}}{1 + \sqrt{\varepsilon}} \right|^2$$
- **Rough Surface (Choudhury et al. 1982):**
  Random micro-roughness reduces coherent specular reflection and increases emissivity:
  $$e_r = 1 + (e_s - 1) e^h$$
  where $h = 4 k_0^2 s^2$ ($s$ is the surface RMS height).
- **Brightness Temperature (Eq. 1.31):**
  $$T_B = \tau_{\text{atm}} [ (1 - e) T_{\text{sky}} + e T_{\text{soil}} ] + (1 - \tau_{\text{atm}}) T_{\text{atm}}$$

---

## 5. Saline Water Dielectric Properties (Behari 2005, Ch. 2)

Dissolved salts dissociate into ions ($\text{Na}^+$, $\text{Cl}^-$) that drift in phase with the low-frequency electric field, producing significant Ohmic conduction:

$$\varepsilon_{\text{saline}} = \varepsilon_{\text{Debye}} + j \frac{\sigma}{\varepsilon_0 \omega}$$

1. **Static Permittivity $\varepsilon_s(S, T)$ (Klein & Swift 1977, Eq. 2.19):**
   $$\varepsilon_s(S) = \varepsilon_{s,0}(T) [1 - 0.2551 S_{\text{frac}} + 5.151\times 10^{-3} S_{\text{frac}}^2 - 6.889\times 10^{-5} S_{\text{frac}}^3]$$
2. **Ionic Conductivity $\sigma(S, T)$ (Stogryn 1971, Eq. 2.21):**
   Empirically parameterized using salinity in parts per thousand ($S_{\text{ppt}}$). Conduction dominates the imaginary loss at L-band (1.4 GHz).

---

## 6. Advanced Dielectric Mixing Theories (Behari 2005, Ch. 6)

Rather than empirical polynomial fitting, physics-based mixing models describe effective medium permittivity $\varepsilon_{\text{eff}}$ from constituent components:

1. **Wang & Schmugge (1980):**
   Separates water into tightly bound water ($m_v \le m_t$) and free water ($m_v > m_t$) where $m_t = 0.49 W_P + 0.165$.
2. **Dobson Semi-Empirical (1985):**
   $$\varepsilon_{\text{eff}}^\alpha = 1 + \frac{\rho_b}{\rho_{ss}} (\varepsilon_{ss}^\alpha - 1) + m_v^\beta (\varepsilon_{fw}^\alpha - 1)$$
   with $\alpha = 0.65$, $\beta \approx 1.09$, soil solids density $\rho_{ss} \approx 2.65\text{ g/cm}^3$, and dry bulk density $\rho_b$.
3. **Rayleigh Mixing (Eq. 6.24):** Electrostatic dipole formulation for dilute spherical inclusions.
4. **Böttcher / Polder–van Santen (Eq. 6.33):** Self-consistent symmetric effective medium equation.
5. **Complex Refractive Index Model (CRIM / Refractive, Eq. 6.12):** Volume-weighted square-root sum.
6. **4-Component Model (Sahu 1998, Eq. 6.39):** Complete separation of solid mineral particles, bound water, free water, and air voids.

---

## 7. Interferometric Coherence & Phase Triplet (De Zan et al. 2014)

### Complex Coherence (Eq. 12–14)
Under the Born volume scattering assumption, the interferogram between two moisture states (wavenumbers $k'_{z1}$ and $k'_{z2}$) is:

$$I(k'_{z1}, k'_{z2}) = \int_0^\infty e^{-2j(k'_{z1} - k'^*_{z2})z} dz = \frac{1}{2j k'_{z1} - 2j k'^*_{z2}}$$

Normalizing by individual powers yields complex coherence $\gamma$:

$$\gamma = \frac{I(k'_{z1}, k'_{z2})}{\sqrt{I(k'_{z1}, k'_{z1}) \cdot I(k'_{z2}, k'_{z2})}}$$

- Magnitude $|\gamma| \le 1$: Decorrelation caused by changing penetration depth and phase mismatch.
- Phase $\angle \gamma$: Interferometric phase bias induced purely by dielectric change.

### Phase Triplet (Closure Phase, Eq. 15)
For three acquisitions $(1, 2, 3)$:

$$\phi_{123} = \phi_{12} + \phi_{23} - \phi_{13} = \angle(\gamma_{12} \gamma_{23} \gamma_{13}^*)$$

While geometric deformation and atmospheric delays cancel exactly ($\Delta \phi_{\text{def}} = 0$), soil moisture dielectric phase variations are non-linear, resulting in $\phi_{123} \ne 0$.

---

## 8. SLC Simulation & Multi-Temporal Closure Phase (Zheng & Fattahi 2026)

### Discretized Depth Sum (Eq. 3)
A single-look radar pixel at time $t$ over depth $D$ divided into $N$ layers is:

$$s(t) = \sum_{n=1}^N \theta_n \exp\bigl(-2j k'_{z,n}(t) z_n\bigr)$$

where $\theta_n \sim \mathcal{CN}(0, \sigma_n)$ is the complex backscatter amplitude.
- **Uniform Profile:** $\sigma_n = 1/N$.
- **Exponential Profile:** $\sigma_n \propto \exp(-z_n / \ell)$, capturing near-surface compaction or root density.

### Closure Phase Time Series
- **Bandwidth-1 (BW-1):** Sequential nearest-neighbor unwrapped phase $\phi_{\text{BW-1}}(t) = \sum_{i=1}^{t-1} \angle(z_{i, i+1})$.
- **Full-Network:** Least-squares phase inversion using all pairs.
- **Closure Phase:** $\phi_{\text{closure}}(t) = \phi_{\text{BW-1}}(t) - \phi_{\text{FullNet}}(t)$.

---

## 9. Soil Moisture Inversion Algorithms

### 1. Coherence Magnitude Root-Finding (De Zan Sec. III-D)
Given known master moisture $m_{v1}$ and observed coherence $|\gamma_{\text{obs}}|$, solve via 1D root-finding (Brent's method):
$$|\gamma_{\text{model}}(m_{v1}, m_{v2})| - |\gamma_{\text{obs}}| = 0$$

### 2. Closure Phase InSAR-SMI Linear Transfer (Zheng Sec. 4)
The closure phase step $\beta$ is linearly proportional to moisture anomaly $\alpha = \Delta m_v$:
$$\beta \approx T_{\text{coeff}} \times \alpha$$
where $T_{\text{coeff}}$ [deg / ($m^3/m^3$)] is the calibrated transfer coefficient. The instantaneous moisture index is:
$$\text{SMI}(t) = \frac{\nabla_t \phi_{\text{closure}}(t)}{T_{\text{coeff}}}$$

### 3. Joint Multi-Objective Inversion
Minimizes weighted squared residual between observed and forward-modeled complex coherence:
$$\mathcal{L}(m_{v2}) = w_m (|\gamma_{\text{obs}}| - |\gamma(m_{v2})|)^2 + w_p (\angle \gamma_{\text{obs}} - \angle \gamma(m_{v2}))^2$$

---

## 10. Comprehensive Symbol & Parameter Glossary

| Symbol | Description | Units | Default / Typical Value |
| :--- | :--- | :--- | :--- |
| $m_v$ | Volumetric soil moisture | $m^3/m^3$ | $0.02 - 0.45$ |
| $\varepsilon = \varepsilon' + j \varepsilon''$ | Complex relative dielectric permittivity | - | $3.5 + 0.1j$ (dry) to $25 + 6j$ (wet) |
| $f$ / `freq_ghz` | Radar frequency | GHz | 1.4 (L), 5.0 / 5.4 (C), 10.0 (X) |
| $\lambda$ | Free-space radar wavelength ($c/f$) | m | 0.214 m (L-band), 0.056 m (C-band) |
| $\theta_{inc}$ | Radar incidence angle | degrees | $30^\circ - 45^\circ$ |
| $k_0$ | Free-space wavenumber ($2\pi / \lambda$) | rad/m | $\approx 29.3\text{ rad/m}$ at L-band |
| $k'_z$ | Vertical complex wavenumber in soil | rad/m | $\text{Im}(k'_z) < 0$ |
| $\delta$ | One-way electromagnetic penetration depth | m | $0.03 - 0.25\text{ m}$ |
| $S, C$ | Sand and Clay mass fractions | - | $S=0.51, C=0.13$ (sandy loam) |
| $\rho_b$ | Soil dry bulk density | $g/cm^3$ | $1.2 - 1.6$ |
| $e_s, e_r$ | Smooth and rough surface emissivity | - | $0.6 - 0.95$ |
| $T_B$ | Radiometric brightness temperature | Kelvin | $200 - 290\text{ K}$ |
| $|\gamma|$ | Complex interferometric coherence | - | $0.0 - 1.0$ |
| $\phi_{123}$ | InSAR phase triplet / closure phase | degrees | $0^\circ - 45^\circ$ |
| $T_{\text{coeff}}$ | Closure phase transfer coefficient | $\text{deg} / (m^3/m^3)$ | $\approx 50 - 150^\circ$ at L-band |
| $M$ / `n_pixels` | Number of looks (pixels) | - | $\ge 100$ |
| $N$ / `n_layers` | Number of vertical depth integration slices | - | $100 - 200$ |
| $D$ / `max_depth`| Maximum numerical integration depth | m | $0.25 - 0.50\text{ m}$ |
| $\ell$ / `scale_depth` | Scatterer exponential decay depth | m | $0.05 - 0.25\text{ m}$ |
