"""Texture-aware Dobson-family formulation from Ulaby textbook companion.

Source: https://ulabyweb.eecs.umich.edu/mrs/codes/Module4_7/plotLibrary.js
Equations 4.66/4.68; high-frequency conductivity branch only (1.4–18 GHz).
This is a specifically identified textbook variant, not a claim to reproduce
every coefficient/version of Dobson (1985). Legacy behari.dobson_semiemp remains.
"""
import numpy as np

def dobson_ulaby(mv, freq_ghz, sand=.51, clay=.13, rho_b=1.3, temperature_c=20.):
    """Return epsilon real + i loss; texture fractions, density g/cm3.

    Dry limit is evaluated analytically. Out-of-domain negative effective loss
    is rejected, never silently clipped. Temperature/density are site inputs,
    not inferred from texture.
    """
    m,s,c,rho,t=np.broadcast_arrays(*[np.asarray(v,float) for v in
                                    (mv,sand,clay,rho_b,temperature_c)])
    if (not np.isfinite(freq_ghz) or not 1.4<=freq_ghz<=18 or
        not all(np.isfinite(v).all() for v in (m,s,c,rho,t)) or
        np.any((s<0)|(c<0)|(s+c>1)|(rho<=0)|(rho>=2.65)|(t<0)|(t>40)|(m<0)|(m>1-rho/2.65))):
        raise ValueError('Unsupported frequency, texture, porosity, density or temperature')
    f=freq_ghz*1e9; alpha=.65
    beta1=1.27-.519*s-.152*c
    beta2=2.06-.928*s-.255*c
    sigma=-1.645+1.939*rho-2.256*s+1.594*c
    static=88.045-.4147*t+6.295e-4*t*t+1.075e-5*t*t*t
    wt=f*(1.1109e-10-3.824e-12*t+6.938e-14*t*t-5.096e-16*t*t*t)
    water_real=4.9+(static-4.9)/(1+wt*wt)
    water_loss=wt*(static-4.9)/(1+wt*wt)
    real=(1+.66*rho+m**beta1*water_real**alpha-m)**(1/alpha)
    # m**beta2 * (water_loss + porosity*sigma/(m*omega*eps0))
    loss=m**beta2*water_loss+m**(beta2-1)*(1-rho/2.65)*sigma/(2*np.pi*8.854e-12*f)
    if np.any(loss<0):
        raise ValueError('Empirical model predicts negative loss for these inputs')
    return real+1j*loss
