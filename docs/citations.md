# NavaSAR References and Citation Guidance

NavaSAR reproduces and builds upon foundational scientific formulations from:

- **De Zan et al. (2014)**, *IEEE Transactions on Geoscience and Remote Sensing* (Born volume scattering, complex coherence, phase triplets).
- **Zheng & Fattahi (2026)**, *Remote Sensing of Environment* (discretized depth sum, closure phase time-series, linear transfer function, and InSAR-SMI retrieval).
- **Behari (2005)**, *Springer* (water relaxation models, 7 soil dielectric mixing models, surface emissivity, brightness temperature, and saline-water electrodynamics).
- **Hallikainen et al. (1985)**, *IEEE Transactions on Geoscience and Remote Sensing* (empirical polynomial dielectric formulations across microwave radar frequencies).
- **Ulaby & Long (2014)**, *Microwave Radar and Radiometric Remote Sensing* (natural-material dielectric models, radar equations, fading, and polarimetry).
- **Ulaby, Moore & Fung (1986)**, *Microwave Remote Sensing, Vol. III* (active microwave sensing of land).

## Primary literature for citation

If you use NavaSAR in academic research, theses, or publications, please cite the primary literature as follows:

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

## Why Ulaby remains important for SAR phase work

Ulaby's texts are especially valuable for the sign conventions, propagation phase terms, phase-center formulations, radar equation structure, and polarization/compatibility conventions that underpin SAR and InSAR interpretation. In other words, the current implementation already relies heavily on the more recent soil-moisture-specific formulations of De Zan and Zheng, but Ulaby is still the best reference for the underlying electromagnetic and system-level phase bookkeeping.

That makes the Ulaby additions most valuable in three layers:

1. **Documentation and conceptual reference**: describe how phase is defined in SAR, interferometry, and polarimetry, with explicit sign conventions and transfer-function language.
2. **API-level phase utilities**: add small helper functions for phase convention checks, propagation-phase bookkeeping, and calibration metadata, without changing the core retrieval formulations.
3. **Research-grade validation**: compare the package's phase and coherence outputs to canonical Ulaby formulations when the objective is system-level consistency rather than new inversion physics.

## Suggested priority for code/package additions

For NavaSAR, the practical priority is:

- **High priority**: maintain the De Zan/Zheng soil-moisture formulations as the core retrieval backbone.
- **Medium priority**: add Ulaby-centered documentation and symbolic phase conventions in the reference material and API docs.
- **Medium priority**: add utility functions for phase-sign conventions, propagation factors, and radar-geometry bookkeeping in the package.
- **Lower priority**: major new algorithms based on Ulaby alone, unless the project explicitly aims to expand beyond soil-moisture retrieval into broader SAR system physics, land backscatter modeling, or polarimetric inversion.

In short: keep the current scientific core centered on De Zan and Zheng, and use Ulaby as the foundation for phase semantics, derivation clarity, and cross-checking, rather than as the primary source for the retrieval model itself.
