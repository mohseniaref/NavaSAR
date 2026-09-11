# SAR and polarimetry foundations

This chapter introduces the ideas used by `navasar.materials`,
`navasar.polarimetry`, and `navasar.land`. It assumes introductory calculus,
complex numbers, probability, and linear algebra.

## Complex permittivity

A microwave field makes charges inside a material move. Some energy is stored
and returned; some becomes heat:

$$\epsilon_r=\epsilon' + j\epsilon''.$$

$\epsilon'$ controls wave speed and refraction, while $\epsilon''$ represents
loss. NavaSAR uses a positive imaginary part for passive loss. References using
$\epsilon'-j\epsilon''$ have chosen the opposite time convention.

Water has a large, frequency-dependent permittivity. Ice, air and dry minerals
have smaller values. This contrast explains radar sensitivity to moisture.

## Natural materials are mixtures

Snow is air plus ice; wet snow also contains water; sea ice contains brine;
soil contains minerals, air and water; vegetation contains dry matter plus free
and bound water. An effective-medium model replaces this microscopic structure
with one bulk permittivity.

There is no universal mixing rule. The answer depends on fractions, inclusion
shapes, orientations and interactions:

- Dilute ellipsoid models suit sparse inclusions.
- Polder–van Santen treats components symmetrically.
- TVB represents inclusions and host as confocal shapes.
- Refractive mixing averages the complex refractive index.

An ellipsoid has factors $L_a+L_b+L_c=1$. A sphere has three factors of $1/3$;
discs and needles respond differently. Shape therefore matters even when
materials and fractions are unchanged.

Pure ice has $\epsilon'\approx3.2$ through much of the microwave band. Air
reduces dry-snow permittivity, whereas a little liquid water strongly increases
wet-snow loss. Sea-ice salt resides in liquid brine pockets whose salinity and
volume change with temperature. Powdering rock introduces air, making bulk
density important. Vegetation contains both fast-relaxing free water and
slow-relaxing bound water.

`vegetation_permittivity` describes plant material.
`canopy_effective_permittivity` mixes that material with air. They represent
different physical scales.

## Scattering matrix and polarization

A polarimetric radar transmits and receives horizontal and vertical fields:

$$S=\begin{bmatrix}S_{HH}&S_{HV}\\S_{VH}&S_{VV}\end{bmatrix}.$$

Each element is complex. For a reciprocal monostatic target,
$S_{HV}=S_{VH}$ under compatible conventions. For transmit state $e_t$ and
receive state $e_r$, polarization synthesis gives

$$s=e_r^H S e_t.$$

Forward-scatter alignment (FSA) and backscatter alignment (BSA) define bases
relative to propagation differently. NavaSAR changes the receive V sign when
converting. Convention errors may reverse cross-polarized phases without
changing power, so convention metadata is essential.

## Radar equation and normalized backscatter

For a monostatic point target,

$$P_r=\frac{P_tG^2\lambda^2\sigma}{(4\pi)^3R^4L}.$$

The $R^4$ term represents spreading on both paths. Point-target radar cross
section $\sigma$ has units of square metres. Distributed terrain uses
dimensionless $\sigma^0$; conversion to RCS requires projected illuminated area.

## Speckle and multilooking

Fields from many unresolved scatterers add coherently. Small phase differences
produce the bright and dark interference called speckle. In the classical
fully-developed model, complex amplitude is circular Gaussian, amplitude is
Rayleigh distributed and single-look intensity is exponential.

An average of $N$ independent intensity looks is Gamma distributed:

$$E[I]=P,\qquad \operatorname{var}(I)=P^2/N.$$

Multilooking reduces speckle at the cost of resolution. A Lee filter uses local
statistics to smooth uniform areas while retaining strong edges; it does not
recover information destroyed by interference.

For repeated complex observations, coherent power is $|E[S]|^2$ and incoherent
power is $E[|S-E[S]|^2]$. Stable specular scattering contributes to the first;
fluctuating rough-surface or volume scattering contributes to the second.

## Covariance, coherency and eigenvalues

For scattering vector $k$,

$$C=E[kk^H].$$

$C$ is Hermitian and positive semidefinite. Diagonal entries are powers;
off-diagonal entries are complex correlations. A Pauli vector reorganizes
channels into useful scattering combinations, and its covariance is the
coherency matrix $T$.

Normalized eigenvalues $p_i$ describe relative orthogonal mechanisms.
Polarimetric entropy and anisotropy are

$$H=-\sum_{i=1}^3p_i\log_3p_i,\qquad
A=\frac{\lambda_2-\lambda_3}{\lambda_2+\lambda_3}.$$

$H=0$ means one mechanism dominates; $H=1$ means equal weights. These numbers
summarize diversity but do not identify land cover without further assumptions.

## Common mistakes

- Mixing dB and linear power inside equations.
- Confusing plant-material and canopy-scale permittivity.
- Confusing gravimetric and volumetric moisture.
- Ignoring dielectric sign or FSA/BSA conventions.
- Treating speckle as additive Gaussian image noise.
- Applying empirical models outside documented ranges.

