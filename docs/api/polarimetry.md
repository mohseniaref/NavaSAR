# Polarimetric radar

`navasar.polarimetry` implements the convention-aware algebra and statistical
tools needed by Chapter 5:

- scattering matrices and Jones polarization states;
- FSA/BSA conversion, radar equation and distributed-target normalization;
- arbitrary transmit/receive polarization synthesis;
- reciprocal Pauli target vectors;
- covariance and Pauli coherency matrices;
- Hermitian eigendecomposition, entropy, and anisotropy;
- correlated circular complex-Gaussian speckle simulation.
- Rayleigh and N-look statistics, coherent/incoherent power and Lee filtering;
- compact span, eigenvalue, entropy and anisotropy parameters.

The matrix channel order is H,V. `synthesize_scattering` evaluates
$e_r^H S e_t$. Covariance matrices use $C=\langle kk^H\rangle$.
