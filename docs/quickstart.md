# Quick start

This guide starts from a new checkout. You do not need an environment named
`code`; that was only the developer's local environment name.

## 1. Create an isolated Python environment

Choose either Conda or standard Python. Do not run both sets of commands.

### Conda

```bash
conda create -n navasar-env python=3.10 -y
conda activate navasar-env
```

### Python `venv`

```bash
python3 -m venv .venv
source .venv/bin/activate
```

On Windows PowerShell: `.venv\Scripts\Activate.ps1`.

## 2. Install NavaSAR

From the directory containing `pyproject.toml`, run:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[notebooks]"
```

Confirm which environment and package Python sees:

```bash
python -c "import sys, navasar; print(sys.executable); print(navasar.__file__)"
```

If an older checkout reports “multiple top-level packages discovered,” update
it. Current NavaSAR explicitly declares `navasar` in `pyproject.toml`.

## 3. Start Jupyter

From the repository root:

```bash
python -m jupyter notebook
```

Using `python -m jupyter` prevents a Jupyter executable from another Python
installation being selected. If a kernel cannot import NavaSAR, register the
active environment and select it in Jupyter:

```bash
python -m ipykernel install --user --name navasar-env --display-name "Python (NavaSAR)"
```

## Notebook learning path

1. Notebooks 01–06: complex SAR signals, coherence, phase and vegetation.
2. Notebooks 07–08: water and soil dielectric behavior.
3. Notebooks 09–11: sensitivity and soil-moisture retrieval.
4. Notebook 12: land backscatter, canopy attenuation and speckle.
5. Notebook 13: natural-material dielectric mixtures.
6. Notebook 14: polarimetric scattering and covariance.

## Minimal SAR and InSAR example

```python
import numpy as np
from navasar.dielectric import hallikainen_dielectric
from navasar.coherence import interferometric_coherence
from navasar.slc import simulate_slc

eps = hallikainen_dielectric(0.25, freq_ghz=1.4, sand=0.51, clay=0.13)
gamma = interferometric_coherence(0.10, 0.25, freq_ghz=1.4)
slc = simulate_slc(np.array([0.10, 0.25]), freq_ghz=1.4)
print(eps, abs(gamma), slc.shape)
```

## Material and polarimetric example

```python
from navasar.materials import pure_ice_permittivity, vegetation_permittivity
from navasar.polarimetry import scattering_matrix, polarization_vector, synthesize_scattering

ice = pure_ice_permittivity(5.4, temperature_c=-10)
leaf = vegetation_permittivity(5.4, gravimetric_moisture=0.7)
S = scattering_matrix(shh=1, shv=0.1j, svh=0.1j, svv=0.4)
h = polarization_vector(orientation_deg=0)
hh = synthesize_scattering(S, transmit=h, receive=h)
```

## Tests

```bash
python -m pip install -e ".[test,notebooks]"
python -m pytest -q
```

