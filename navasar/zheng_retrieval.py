"""Zheng & Fattahi (2026), Eqs. 10--11 and Section 6.1.

Event fitting is an independent numerical implementation of the published
equations, not the authors' software. BIC, optimization bounds and event-search
rules below are implementation choices. Outputs are relative, not VWC.
"""
import numpy as np
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from scipy.signal import find_peaks


def normalize_smi(values, minimum_range=1e-10):
    """Min/max normalization, not clipping. Constant/missing data have no index."""
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)) or np.ptp(values) <= minimum_range:
        return np.full_like(values, np.nan)
    return (values-values.min())/np.ptp(values)


def event_response(time, beta, omega, onset):
    """Eq. 10 closure and Eq. 11 unnormalized moisture, A=lambda=1.

    Units: beta is degrees, omega and onset have the units of time; moisture
    output is therefore an index precursor (degrees/time), NOT m3/m3.
    """
    t = np.asarray(time, dtype=float)[:, None]
    beta, omega, onset = map(np.atleast_1d, (beta, omega, onset))
    if np.any(omega <= 0):
        raise ValueError('Event time constants must be positive')
    elapsed = np.maximum(t-onset, 0)
    active = t >= onset
    closure = np.sum(beta*(-np.expm1(-elapsed/omega))*active, axis=1)
    moisture = np.sum((beta/omega)*np.exp(-elapsed/omega)*active, axis=1)
    return closure, moisture


def retrieve_smi(time, closure_deg, max_events=6, median_window=3,
                 minimum_range_deg=1e-6, loss='linear'):
    """Fit signed saturating steps and return both paper-based SMI estimators.

    Uses actual acquisition times; never interpolates gaps or fills NaNs.
    The closure gradient proposes onset neighborhoods, BIC selects event count
    including a constant/no-event null. Local onset and decay are optimized.
    A=lambda=1 follows Section 4. Median filtering follows Section 5.
    Robust soft_l1 is an optional extension, not claimed as the paper's solver.
    """
    t, y = np.asarray(time, float), np.asarray(closure_deg, float)
    if t.ndim != 1 or y.shape != t.shape or len(t) < 8:
        raise ValueError('Need matching one-dimensional arrays with >=8 dates')
    if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError('Time must be finite and strictly increasing')
    if median_window < 1 or median_window % 2 != 1 or median_window > len(t):
        raise ValueError('Median window must be odd and <= number of dates')
    if loss not in ('linear', 'soft_l1') or max_events < 1:
        raise ValueError('Invalid loss or max_events')
    nan = np.full_like(y, np.nan)
    if not np.all(np.isfinite(y)):
        return dict(status='missing_closure', event_smi=nan, gradient_smi=nan)
    smooth = median_filter(y, size=median_window, mode='nearest')
    gradient = np.gradient(smooth, t)
    result = dict(time=t, observed=y, smoothed=smooth, gradient=gradient,
                  gradient_smi=normalize_smi(gradient), event_smi=nan,
                  loss=loss, A=1., temporal_scale_lambda=1.)
    if np.ptp(smooth) <= minimum_range_deg:
        result.update(status='no_resolved_closure', gradient_smi=nan)
        return result
    n = len(t)
    dt = float(np.median(np.diff(t)))
    span = float(t[-1]-t[0])
    scale = max(float(np.ptp(smooth)), minimum_range_deg)
    # Robust scale for optional M-estimation; all candidates use the same scale.
    noise = max(float(np.median(np.abs(y-smooth))), scale*.02, minimum_range_deg)
    score = np.abs(gradient)
    peaks = list(find_peaks(score, distance=2)[0])
    # Boundary events remain candidates, rather than manufacturing a baseline.
    peaks += [1, n-2]
    ranked = sorted(set(peaks), key=lambda k: score[k], reverse=True)
    ranked = [k for k in ranked if score[k] > 0]
    null = np.full(n, np.mean(smooth))
    def bic(pred, parameters):
        mse = max(np.mean((pred-smooth)**2), minimum_range_deg**2)
        return float(n*np.log(mse)+parameters*np.log(n))
    best = dict(bic=bic(null, 1), predicted=null, events=[], count=0)
    candidates = [dict(count=0, bic=best['bic'])]
    for count in range(1, min(max_events, len(ranked), (n-2)//3)+1):
        selected = sorted(ranked[:count])
        lower, upper, initial = [smooth[0]-2*scale], [smooth[0]+2*scale], [smooth[0]]
        for k in selected:
            lo, hi = t[max(0,k-2)], t[min(n-2,k+1)]
            if hi <= lo:
                hi = lo+dt*.5
            lower += [-10*scale, max(dt*.25, 1e-6), lo]
            upper += [10*scale, span*2, hi]
            initial += [np.sign(gradient[k])*scale/count, max(dt,span/12),
                        np.clip(t[max(k-1,0)], lo+1e-8, hi-1e-8)]
        def predict(x):
            pars = np.asarray(x[1:]).reshape(-1,3)
            return x[0]+event_response(t, *pars.T)[0]
        fit = least_squares(lambda x: (predict(x)-smooth)/noise, initial,
                            bounds=(lower,upper), loss=loss, max_nfev=600)
        pred = predict(fit.x)
        cost = bic(pred, len(fit.x))
        candidates.append(dict(count=count, bic=cost, converged=bool(fit.success)))
        if fit.success and cost < best['bic']:
            best = dict(bic=cost, predicted=pred, count=count,
                        events=fit.x[1:].reshape(-1,3), offset=float(fit.x[0]))
    result.update(closure_fit=best['predicted'], candidate_scores=candidates,
                  events=np.asarray(best['events']).reshape(-1,3),
                  event_count=best['count'], bic=best['bic'],
                  closure_rmse_deg=float(np.sqrt(np.mean((best['predicted']-smooth)**2))))
    if best['count']:
        precursor = event_response(t, *np.asarray(best['events']).T)[1]
        result.update(status='relative_smi', event_precursor=precursor,
                      event_smi=normalize_smi(precursor),
                      omega_lower_bound_count=int(np.sum(np.asarray(best['events'])[:,1]<=dt*.25*1.01)),
                      max_event_limit_reached=bool(best['count']==max_events))
    else:
        result.update(status='no_events_preferred')
    return result
