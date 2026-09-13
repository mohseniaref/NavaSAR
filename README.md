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

For the foundational references, citation guidance, and a research-oriented overview of the scientific lineage behind NavaSAR, see [docs/citations.md](docs/citations.md).

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

See [docs/repository_structure.md](docs/repository_structure.md) for the full directory tree and a description of every package module, notebook, and documentation file.

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

## 📜 License

NavaSAR is distributed under the open-source **[MIT License](LICENSE)**. You are free to use, modify, and distribute the code for academic, research, and commercial purposes with attribution.
