<div align="center">

# NavaSAR

### Physics-Based SAR and InSAR Modeling, Simulation & Retrieval

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Release: alpha](https://img.shields.io/badge/release-0.1.0a1-orange.svg)](https://github.com/mohseniaref/NavaSAR/releases)
[![Tests](https://github.com/mohseniaref/NavaSAR/actions/workflows/tests.yml/badge.svg)](https://github.com/mohseniaref/NavaSAR/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22700038.svg)](https://doi.org/10.5281/zenodo.22700038)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-reference-green.svg)](docs/index.md)

</div>

---

## 🌍 Overview

**NavaSAR** is an open-source Python physics library for research and education in radar remote sensing. It covers natural-material dielectric properties, polarimetric scattering and statistics, detected land imagery, complex SAR simulation, interferometric coherence, closure phase, and soil-moisture retrieval.

NavaSAR reproduces and builds upon foundational scientific formulations from:
- **De Zan et al. (2014)**, *IEEE Transactions on Geoscience and Remote Sensing* (Born volume scattering, complex coherence, phase triplets).
- **Zheng & Fattahi (2026)**, *Remote Sensing of Environment* (Discretized depth sum, closure phase time-series, linear transfer function, and InSAR-SMI retrieval).
- **Behari (2005)**, *Springer* (Water relaxation models, 7 soil dielectric mixing models, surface emissivity, brightness temperature, and saline water electrodynamics).
- **Hallikainen et al. (1985)**, *IEEE Transactions on Geoscience and Remote Sensing* (Empirical polynomial dielectric formulations across microwave radar frequencies).
- **Ulaby & Long (2014)**, *Microwave Radar and Radiometric Remote Sensing* (natural-material dielectric models, radar equations, fading, and polarimetry).
- **Ulaby, Moore & Fung (1986)**, *Microwave Remote Sensing, Vol. III* (active microwave sensing of land).

---

## 🔬 Core Physics at a Glance

1. **Dielectric Permittivity ($\varepsilon = \varepsilon' + j\varepsilon''$):**  
   Water has high relative permittivity ($\varepsilon \approx 80$) relative to dry soil minerals ($\varepsilon \approx 3 - 4$). Increasing moisture $m_v$ increases dielectric storage ($\varepsilon'$) and Ohmic loss ($\varepsilon''$).
2. **Vertical Wavenumber ($k'_z$):**  
   Refraction into the lossy soil medium is governed by $k'_z = \sqrt{k_0^2 \varepsilon_r - k_x^2}$, with the physical branch $\text{Im}(k'_z) < 0$ enforcing downward attenuation.
3. **Volume Scattering & Complex Coherence ($\gamma$):**  
   Integrating radar echoes over depth creates interferometric decorrelation ($|\gamma| < 1$) and phase delay ($\angle \gamma$).
4. **Non-Zero Closure Phase ($\phi_{123} \ne 0$):**  
   Unlike linear ground displacement and atmospheric delays that cancel across closed interferometric loops, non-linear dielectric shifts create persistent non-closure signatures that serve as a direct diagnostic tracer for soil moisture.
5. **Soil Moisture Retrieval (InSAR-SMI):**  
   Moisture anomalies are retrieved quantitatively via 1D coherence root-finding and the temporal gradient of InSAR closure phase time-series ($\text{SMI}(t) = \nabla_t \phi_{\text{closure}}(t) / T_{\text{coeff}}$).

---

## 📁 Repository Structure

```
NavaSAR/
├── navasar/                            # Core Scientific Python Package
│   ├── __init__.py                     # Package exports
│   ├── dielectric.py                   # Hallikainen (1985) permittivity & kz
│   ├── coherence.py                    # De Zan (2014) coherence, triplets, Fresnel phase
│   ├── slc.py                          # Zheng & Fattahi (2026) SLC depth sum simulation
│   ├── closure.py                      # Multi-temporal closure phase (BW-1 vs Full-Network)
│   ├── behari.py                       # Behari (2005) water, mixing, emissivity & saline models
│   ├── penetration.py                  # One-way and two-way radar penetration depth models
│   ├── inversion.py                    # Coherence root-finding & InSAR-SMI retrieval
│   ├── land.py                         # Ulaby Ch. 21 land backscatter & detected-image statistics
│   ├── materials.py                    # Long & Ulaby Ch. 4 material dielectric models
│   └── polarimetry.py                  # Long & Ulaby Ch. 5 polarimetric radar tools
├── notebooks/                          # 14 Curated Interactive Jupyter Notebooks
│   ├── 01_complex_slc_simulation.ipynb           # Dielectric model & kz (De Zan Fig. 2)
│   ├── 02_soil_moisture_effect.ipynb             # Coherence magnitude & phase grids (Fig. 3, 4)
│   ├── 03_deformation_and_moisture.ipynb         # Fresnel transmission & phase triplets (Fig. 5, 8)
│   ├── 04_multilooking_and_fading.ipynb          # Multi-looking & AGRISAR 2006 validation (Fig. 6)
│   ├── 05_closure_phase.ipynb                    # Closure phase time-series (Zheng Fig. 4, 5, 6)
│   ├── 06_vegetation_and_mixed_materials.ipynb   # Polarimetric phase & vegetation decorrelation
│   ├── 07_water_dielectric_models.ipynb          # Debye, Cole-Cole, and Saline water electrodynamics
│   ├── 08_soil_dielectric_mixing_models.ipynb    # 7 Soil mixing models, emissivity & penetration
│   ├── 09_zheng_simulation_retrieval.ipynb       # Anomaly transfer function (Fig. 8) & multi-band
│   ├── 10_zheng_figs7_9_10_sensitivity.ipynb     # End-to-end pipeline (Fig. 7), depth D & multilook noise
│   ├── 11_inversion_and_retrieval.ipynb          # Quantitative retrieval: Brent root-finding & InSAR-SMI
│   ├── 12_land_backscatter_simulation.ipynb      # Canopy, row direction & multilook land imagery
│   ├── 13_material_dielectric_models.ipynb       # Ellipsoid mixtures and dry snow
│   └── 14_polarimetric_sar_and_speckle.ipynb     # Polarization synthesis and covariance
├── docs/                               # Comprehensive Documentation & Undergrad Primer
│   ├── index.md                        # Documentation overview
│   ├── theory.md                       # Complete mathematical derivation & symbol glossary
│   ├── quickstart.md                   # Installation & usage walkthrough
│   ├── figures.md                      # Index of all 24+ reproduced scientific figures
│   └── api/                            # Detailed API references
├── examples/                           # Auto-generated high-resolution figures
├── pyproject.toml                      # Modern PEP 517/621 package build configuration
├── CITATION.cff                        # Machine-readable academic citation metadata
├── .zenodo.json                        # Zenodo repository archiving configuration
├── LICENSE                             # MIT Open-Source License
└── README.md
```

---

## 🚀 Quick Start

### 1. Installation

Clone the repository and install in development mode:

```bash
git clone https://github.com/mohseniaref/NavaSAR.git
cd NavaSAR

# Install into your active environment
python -m pip install -e ".[notebooks]"
```

### 2. Interactive Notebooks

Launch Jupyter from the repository root and explore notebooks `01` through `14`:

```bash
python -m jupyter notebook
```

### 3. Programmatic Usage

```python
import numpy as np
from navasar.dielectric import hallikainen_dielectric, kz_soil
from navasar.coherence import interferometric_coherence, phase_triplet
from navasar.slc import simulate_slc
from navasar.closure import closure_phase_timeseries
from navasar.inversion import calibrate_closure_transfer, insar_smi_timeseries

# 1. Forward dielectric permittivity (Hallikainen 1985)
eps = hallikainen_dielectric(mv=0.22, freq_ghz=1.4, sand=0.51, clay=0.13)

# 2. Interferometric coherence between two moisture states (De Zan 2014)
gamma = interferometric_coherence(mv1=0.08, mv2=0.22, freq_ghz=1.4)
print(f"|gamma| = {abs(gamma):.3f}, phase = {np.rad2deg(np.angle(gamma)):.1f} deg")

# 3. Simulate multi-temporal SLC stack with rainfall event (Zheng & Fattahi 2026)
mv_ts = np.array([0.05]*10 + [0.20]*4 + [0.05]*10)
slc_stack = simulate_slc(mv_ts, freq_ghz=1.4, sigma_profile='exponential', scale_depth=0.15)

# 4. InSAR closure phase time-series
cp_ts, bw1_ts, fullnet_ts = closure_phase_timeseries(slc_stack)

# 5. Retrieve Soil Moisture Index (InSAR-SMI)
_, _, T_coeff = calibrate_closure_transfer(mv_base=0.05, freq_ghz=1.4)
smi = insar_smi_timeseries(cp_ts, transfer_coeff=T_coeff)
```

---

## 📖 Key References & Citations

If you utilize NavaSAR in academic research, theses, or publications, please cite the primary literature:

- **Ulaby, F. T., & Long, D. G. (2014).**  
  *Microwave Radar and Radiometric Remote Sensing.*  
  University of Michigan Press, ISBN 978-0-472-11935-6.

- **Ulaby, F. T., Moore, R. K., & Fung, A. K. (1986).**  
  *Microwave Remote Sensing: Active and Passive, Volume III—From Theory to Applications.*  
  Artech House.

- **De Zan, F., Parizzi, A., Prats-Iraola, P., & López-Dekker, P. (2014).**  
  *A SAR interferometric model for soil moisture.*  
  **IEEE Transactions on Geoscience and Remote Sensing**, 52(1), 418–425.  
  DOI: [10.1109/TGRS.2013.2241069](https://doi.org/10.1109/TGRS.2013.2241069)

- **Zheng, Y., & Fattahi, H. (2026).**  
  *Modeling, prediction, and retrieval of surface soil moisture from InSAR closure phase.*  
  **Remote Sensing of Environment**, 333, 115104.  
  DOI: [10.1016/j.rse.2025.115104](https://doi.org/10.1016/j.rse.2025.115104)

- **Behari, J. (2005).**  
  *Microwave Dielectric Behavior of Wet Soils.*  
  **Springer Science & Business Media**, ISBN: 1-4020-3271-4.  
  DOI: [10.1007/1-4020-3282-9](https://doi.org/10.1007/1-4020-3282-9)

- **Hallikainen, M. T., Ulaby, F. T., Dobson, M. C., El-Rayes, M. A., & Wu, L.-K. (1985).**  
  *Microwave dielectric behavior of wet soil – Part 1: Empirical models and experimental observations.*  
  **IEEE Transactions on Geoscience and Remote Sensing**, GE-23(1), 25–34.  
  DOI: [10.1109/TGRS.1985.289497](https://doi.org/10.1109/TGRS.1985.289497)

---

## 📜 License

NavaSAR is distributed under the open-source **[MIT License](LICENSE)**. You are free to use, modify, and distribute the code for academic, research, and commercial purposes with attribution.
