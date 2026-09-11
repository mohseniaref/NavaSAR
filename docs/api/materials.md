# Material dielectric models

`navasar.materials` collects reusable Chapter 4 dielectric building blocks.
It reuses the existing water implementations and adds ice, snow, sea ice,
rock, vegetation, and general effective-medium relations.

Model groups include water and ice, dry and wet snow, sea ice, powdered rock,
vegetation material and canopy-scale mixtures. General mixture tools include
ellipsoid depolarization factors, dilute mixing, Polder–van Santen, TVB, and
refractive-index formulations.

Complex values use the NavaSAR loss convention $\epsilon=\epsilon'+j\epsilon''$.
Mixture fractions are volume fractions and must sum to one where applicable.
The compact sea-ice and wet-snow functions remain effective-medium models and
do not replace ice-type-specific field calibration.
