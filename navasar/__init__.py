"""
navasar: SAR interferometric soil moisture physics library.

Implements models from:
  De Zan et al. (2014), IEEE TGRS 52(1), 418-425
  Zheng & Fattahi (2026), RSE 333, 115104
  Hallikainen et al. (1985), IEEE TGRS GE-23(1), 25-34
  Behari, J. (2005). Microwave Dielectric Behavior of Wet Soils. Springer.

Public API
----------
Dielectric models (Hallikainen 1985 + De Zan kz):
    hallikainen_dielectric, kz_soil

Interferometric observables (De Zan 2014, Eq. 12-15):
    interferometric_coherence, phase_triplet, fresnel_phase

SLC simulation (Zheng & Fattahi 2026, Eq. 3-4):
    simulate_slc, multilook_interferogram
    sigma_profile : 'uniform' (default) or 'exponential' scatterer profile

Closure phase time-series (Zheng & Fattahi 2026):
    closure_phase_timeseries

Behari (2005) water and soil dielectric models:
    navasar.behari.*              -- includes saline_water_dielectric

Penetration depth (De Zan Eq. 5 + Behari Eq. 1.23):
    navasar.penetration.*

Moisture retrieval (De Zan Sec. III-D + Zheng Sec. 4):
    navasar.inversion.*

Land backscatter and detected-image simulation (Ulaby et al., Ch. 21):
    navasar.land.*
"""
__version__ = "0.1.0a1"

# --- De Zan / Hallikainen core ---
from .dielectric import hallikainen_dielectric, kz_soil

# --- Interferometric observables ---
from .coherence import interferometric_coherence, phase_triplet, fresnel_phase

# --- SLC simulation ---
from .slc import simulate_slc, multilook_interferogram

# --- Closure phase ---
from .closure import closure_phase_timeseries

# --- Behari (2005) models as a sub-namespace ---
from . import behari

# --- Penetration depth ---
from . import penetration

# --- Moisture retrieval / inversion ---
from . import inversion

# --- Land backscatter, vegetation, periodic rows, and image statistics ---
from . import land

# --- Long & Ulaby Chapters 4 and 5 ---
from . import materials
from . import polarimetry
