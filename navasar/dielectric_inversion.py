"""Experimental dielectric-model closure inversion, separate from Zheng SMI.

Uses source-reviewed Hallikainen and Ulaby/Dobson-family implementations;
Wang–Schmugge retains the local implementation and explicit assumed parameters.
These implementations have documented limitations: explicit opt-in is required.
Uniform semi-infinite volume scattering; fixed texture, density and temperature;
one externally specified moisture reference; bounded, multistart optimization.
See docs/soil_moisture_retrieval.md before scientific interpretation.
"""
import numpy as np
from scipy.optimize import least_squares
from .dielectric import hallikainen_dielectric
from .behari import wang_schmugge
from .dobson_reviewed import dobson_ulaby

MODEL_NOTES = {
    'hallikainen': 'Reviewed 2026-10-15 against MIMICS Table E.1: percent-weight texture; 1.4–18 GHz, adjacent interpolation.',
    'dobson': 'Ulaby Module 4.7 high-frequency Dobson-family variant; texture-aware, assumed density/temperature.',
    'wang-schmugge': 'Local empirical implementation with assumed density, temperature and loss parameters.',
}


def dielectric_covariance(moisture, *, model, frequency_ghz, sand, clay,
                          incidence_deg=40., bulk_density=1.3, temperature_c=20.,
                          allow_experimental=False):
    """Normalized complex covariance from a named local dielectric model.

    Moisture/sand/clay are fractions; frequency GHz, density g/cm3, angle degrees.
    Returns C_ij = integral exp(-2ik_i z) exp(+2ik_j* z) dz, normalized.
    Uses NavaSAR's negative-imaginary wavenumber root/sign convention.
    This is the uniform semi-infinite model, not the paper's finite-layer model.
    """
    if not allow_experimental:
        raise ValueError('Read MODEL_NOTES and explicitly set allow_experimental=True')
    if model not in MODEL_NOTES:
        raise ValueError(f'model must be one of {tuple(MODEL_NOTES)}')
    m=np.asarray(moisture,float)
    if m.ndim!=1 or len(m)<3 or not np.all(np.isfinite(m)) or np.any((m<0)|(m>.6)):
        raise ValueError('Need >=3 finite moisture fractions between 0 and 0.6')
    if not (np.isfinite(sand) and np.isfinite(clay) and 0<=sand<=1 and 0<=clay<=1 and sand+clay<=1):
        raise ValueError('Sand/clay must be scalar fractions summing to <=1')
    if not (np.isfinite(frequency_ghz) and frequency_ghz>0 and 0<=incidence_deg<90
            and 0<bulk_density<2.65 and np.isfinite(temperature_c)):
        raise ValueError('Invalid frequency, incidence, density or temperature')
    if model=='hallikainen':
        eps=hallikainen_dielectric(m,frequency_ghz,sand,clay)
    elif model=='dobson':
        eps=dobson_ulaby(m,frequency_ghz,sand,clay,bulk_density,temperature_c)
    else:
        fun=wang_schmugge
        eps=fun(m,frequency_ghz,sand_pct=100*sand,clay_pct=100*clay,
                rho_b=bulk_density,T_celsius=temperature_c)
    k=2*np.pi*frequency_ghz/.299792458*np.sqrt(eps-np.sin(np.deg2rad(incidence_deg))**2)
    k=np.where(k.imag>0,-k,k)
    if not np.all(np.isfinite(k)) or np.any(k.imag>=0):
        raise ValueError('Non-attenuating or non-finite dielectric response')
    cov=1/(2j*(k[:,None]-k[None,:].conj()))
    power=cov.diagonal().real
    return cov/np.sqrt(power[:,None]*power[None,:])


def triangle_closure(covariance, triangles):
    """Wrapped closure in radians, orientation i->j->k->i."""
    i,j,k=np.asarray(triangles,int).T
    return np.angle(covariance[i,j]*covariance[j,k]*covariance[k,i])


def invert_closure_moisture(covariance, time_days, *, reference_moisture,
                            model, frequency_ghz, sand, clay, triangles=None,
                            bounds=(.005,.45), smoothness=.1, min_coherence=.2,
                            starts=3, seed=42, max_nfev=150,
                            allow_experimental=False, **forward_options):
    """Conditional dielectric inversion of wrapped triangle closures.

    reference_moisture fixes acquisition 0, and must be independently justified.
    Objective: weighted circular phase residual /0.1 rad plus smoothness-weighted
    moisture increments /0.1 fraction, scaled by relative time spacing.
    The residual scales and smoothness are explicit research choices, not a
    likelihood calibrated to real observation errors. Triplets are correlated.
    Returned moisture is conditional on the reference, model and regularization;
    bounds or optimizer convergence do not establish absolute accuracy.
    """
    cov=np.asarray(covariance); t=np.asarray(time_days,float); n=len(t)
    if cov.shape!=(n,n) or n<4 or not np.iscomplexobj(cov) or not np.all(np.isfinite(cov)):
        raise ValueError('Need finite complex covariance matching >=4 acquisition times')
    if not np.allclose(cov,cov.conj().T) or np.any(cov.diagonal().real<=0):
        raise ValueError('Expected Hermitian covariance with positive powers')
    if np.any(np.diff(t)<=0) or not np.all(np.isfinite(t)):
        raise ValueError('Time must be finite and strictly increasing')
    if not (0<=bounds[0]<reference_moisture<bounds[1]<=.6) or starts<2 or smoothness<0 or not 0<=min_coherence<=1:
        raise ValueError('Invalid bounds/reference, starts, smoothness or coherence threshold')
    power=cov.diagonal().real
    cov=cov/np.sqrt(power[:,None]*power[None,:])
    if triangles is None:
        triangles=np.array([(i,i+1,k) for i in range(n-2) for k in range(i+2,min(n,i+7))])
    triangles=np.asarray(triangles)
    if triangles.ndim!=2 or triangles.shape[1]!=3 or not np.issubdtype(triangles.dtype,np.integer):
        raise ValueError('triangles must have integer shape (Q,3)')
    if np.any(triangles<0) or np.any(triangles>=n) or np.any(np.diff(np.sort(triangles,axis=1),axis=1)==0):
        raise ValueError('Invalid triangle indices')
    i,j,k=triangles.T
    weight=np.minimum.reduce([abs(cov[i,j]),abs(cov[j,k]),abs(cov[k,i])])
    valid=weight>=min_coherence
    triangles=triangles[valid]; weight=weight[valid]
    if len(triangles)<n-1:
        raise ValueError('Too few coherent triangle observations; no inversion performed')
    obs=triangle_closure(cov,triangles)
    kwargs=dict(model=model,frequency_ghz=frequency_ghz,sand=sand,clay=clay,
                allow_experimental=allow_experimental,**forward_options)
    # Fail before optimization when the model/physical inputs are unsupported.
    dielectric_covariance(np.full(n,reference_moisture),**kwargs)
    interval=np.sqrt(np.diff(t)/np.median(np.diff(t)))
    def residual(x):
        mv=np.r_[reference_moisture,x]
        pred=triangle_closure(dielectric_covariance(mv,**kwargs),triangles)
        circular=np.angle(np.exp(1j*(pred-obs)))
        return np.r_[np.sqrt(weight)*circular/.1,smoothness*np.diff(mv)/(.1*interval)]
    rng=np.random.default_rng(seed); fits=[]
    for _ in range(starts):
        x=np.clip(reference_moisture+rng.normal(0,.03,n-1),bounds[0]+1e-6,bounds[1]-1e-6)
        fits.append(least_squares(residual,x,bounds=bounds,max_nfev=max_nfev))
    best=min(fits,key=lambda f:np.dot(f.fun,f.fun))
    moisture=np.r_[reference_moisture,best.x]
    errors=np.angle(np.exp(1j*(triangle_closure(dielectric_covariance(moisture,**kwargs),triangles)-obs)))
    return dict(moisture=moisture,status='conditional_fit' if best.success else 'iteration_limit',
                model=model,model_note=MODEL_NOTES[model],reference_moisture=reference_moisture,
                model_parameters=kwargs,bounds=bounds,smoothness=smoothness,
                closure_rmse_deg=float(np.rad2deg(np.sqrt(np.mean(errors**2)))),
                triangles=triangles,triangle_count=len(triangles),
                multistart_max_spread=float(np.max(np.ptp(np.array([f.x for f in fits]),axis=0))),
                data_jacobian_rank=int(np.linalg.matrix_rank(best.jac[:len(triangles)])),
                unknowns=n-1,bound_fraction=float(np.mean(np.isclose(best.x,bounds[0],atol=1e-5)|np.isclose(best.x,bounds[1],atol=1e-5))),
                absolute_moisture_validated=False)
