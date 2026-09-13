# NavaSAR Documentation

Physics-based SAR interferometric soil moisture modeling, simulation, and retrieval library.

---

## Documentation Sections

| Page | Description |
|------|-------------|
| [quickstart.md](quickstart.md) | Installation, environment setup, and programmatic workflow |
| [citations.md](citations.md) | Scientific lineage, primary literature, and citation guidance |
| [repository_structure.md](repository_structure.md) | Full directory tree and description of every module, notebook, and doc file |
| [sar_foundations.md](sar_foundations.md) | Undergraduate introduction to material mixtures, polarimetry, radar statistics and conventions |
| [theory.md](theory.md) | Comprehensive undergraduate primer, mathematical derivations & symbol glossary |
| [api/dielectric.md](api/dielectric.md) | `navasar.dielectric` — Hallikainen (1985) permittivity & vertical wavenumber $k'_z$ |
| [api/coherence.md](api/coherence.md) | `navasar.coherence` — De Zan (2014) volume coherence, phase triplets, Fresnel transmission |
| [api/slc.md](api/slc.md) | `navasar.slc` — SLC coherent depth integration & multilooking |
| [api/closure.md](api/closure.md) | `navasar.closure` — Bandwidth-1 vs. Full-Network closure phase time-series |
| [api/behari.md](api/behari.md) | `navasar.behari` — Behari (2005) water, mixing, emissivity & saline models |
| [api/penetration.md](api/penetration.md) | `navasar.penetration` — One-way & two-way radar penetration depth formulations |
| [api/inversion.md](api/inversion.md) | `navasar.inversion` — Coherence root-finding & InSAR-SMI retrieval algorithms |
| [api/land.md](api/land.md) | `navasar.land` — Chapter 21 canopy, row-structure, speckle & detected land-image simulation |
| [api/materials.md](api/materials.md) | `navasar.materials` — Long & Ulaby Ch. 4 dielectric mixtures, ellipsoids & snow |
| [api/polarimetry.md](api/polarimetry.md) | `navasar.polarimetry` — Long & Ulaby Ch. 5 matrices, synthesis, covariance & speckle |
| [figures.md](figures.md) | Catalogue of reproduced figures across literature |

---

## Curated Notebook Tutorials

| Notebook | Focus | Primary Literature Reference |
| :--- | :--- | :--- |
| `01_complex_slc_simulation.ipynb` | Complex dielectric constant & vertical wavevector | De Zan et al. (2014) Fig. 2 |
| `02_soil_moisture_effect.ipynb` | Coherence magnitude & phase grids | De Zan et al. (2014) Figs. 3, 4 |
| `03_deformation_and_moisture.ipynb` | Boundary Fresnel transmission & non-zero phase triplets | De Zan et al. (2014) Figs. 5, 8 |
| `04_multilooking_and_fading.ipynb` | Multilooking speckle reduction & AGRISAR 2006 validation | De Zan et al. (2014) Fig. 6 |
| `05_closure_phase.ipynb` | Closure phase time-series under depth & texture scenarios | Zheng & Fattahi (2026) Figs. 4, 5, 6 |
| `06_vegetation_and_mixed_materials.ipynb` | Polarimetric HH-VV phase & vegetation decorrelation | De Zan et al. (2014) Figs. 9, 10, 11 |
| `07_water_dielectric_models.ipynb` | Debye, Cole-Cole, and Saline water electrodynamics | Behari (2005) Ch. 2 |
| `08_soil_dielectric_mixing_models.ipynb` | 7 Soil dielectric mixing models, emissivity & penetration | Behari (2005) Ch. 1 & 6 |
| `09_zheng_simulation_retrieval.ipynb` | Linear anomaly transfer function & multi-frequency comparison | Zheng & Fattahi (2026) Fig. 8 |
| `10_zheng_figs7_9_10_sensitivity.ipynb` | End-to-end pipeline (Fig. 7), depth $D$ & multilook noise $1/\sqrt{M}$ | Zheng & Fattahi (2026) Figs. 7, 9, 10 |
| `11_inversion_and_retrieval.ipynb` | Quantitative moisture retrieval via coherence & InSAR-SMI | De Zan (2014) & Zheng (2026) |
| `12_land_backscatter_simulation.ipynb` | Canopy attenuation, agricultural rows & N-look detected imagery | Ulaby, Moore & Fung (1986), Ch. 21 |
| `13_material_dielectric_models.ipynb` | Ellipsoidal mixing, Polder–van Santen, generalized power-law mixing, solid rock, dry/wet (Dobson–Peplinski) soil & dry snow | Long & Ulaby, Ch. 4 |
| `14_polarimetric_sar_and_speckle.ipynb` | Polarization synthesis, coherency and correlated speckle | Long & Ulaby, Ch. 5 |
