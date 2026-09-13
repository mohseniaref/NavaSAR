# Repository Structure

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
│   ├── citations.md                    # Scientific lineage, primary literature & citation guidance
│   ├── repository_structure.md         # This file
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
