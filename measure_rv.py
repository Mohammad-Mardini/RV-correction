#!/usr/bin/env python3
# HASHEM: standalone radial-velocity measurement with star-vs-template plots.
#
# Copyright (c) 2026 Mohammad K. H. Mardini
# Licensed under the MIT License; see LICENSE.
# ─────────────────────────────────────────────────────────────────────────────

import argparse
import csv
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, Angle, FK5
import astropy.units as u
from scipy import interpolate
from scipy.interpolate import interp1d
import numpy.polynomial.chebyshev as _cheb
import numpy.polynomial.legendre as _leg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# Settings used by "python measure_rv.py" (relative paths: from the folder you
# run it in).  Command-line arguments override them.
# ─────────────────────────────────────────────────────────────────────────────
SPECTRA  = "*_multi.fits"       # spectra to measure (a glob pattern or a list)
TEMPLATE = "/Users/mohammad/Library/CloudStorage/Dropbox/ref_spectra/hd122563.fits"  # rest-frame RV template
REGION   = "Ca II IR"           # "Ca II IR", "Halpha", "Mg b", "Hbeta", "Ca H&K" or "auto"
OUT_CSV  = "rv_results.csv"     # results table
PLOT_DIR = "rv_plots"           # one figure per spectrum
PLOTS    = True                 # False: measure only

_TEMPLATE = {"path": None}       # template file (used by the joint-CCF fallback)
REGION_ALIASES = {"halpha": "H\u03b1", "hbeta": "H\u03b2"}
CSV_NAMES = {"Ca II IR": "CaIIIR", "H\u03b1": "Halpha", "Mg b": "Mgb",
             "H\u03b2": "Hbeta", "Ca H&K": "CaHK"}


# ─────────────────────────────────────────────────────────────────────────────
# Measurement: copied unchanged from coadd_gui.py
# ─────────────────────────────────────────────────────────────────────────────

CCF_REGIONS   = [(8450,8750,"Ca II IR"),(6510,6610,"Hα"),(5100,5200,"Mg b"),
                 (4810,4910,"Hβ"),(4290,4390,"Ca H&K")]


OBSERVATORY   = EarthLocation(lat=-29.0146*u.deg,lon=-70.6926*u.deg,height=2380*u.m)


C_KMS         = 2.99792458e5


def load_template(path):
    with fits.open(path) as h:
        hdr=h[0].header; flux=h[0].data.astype(np.float64)
        wav=hdr["CRVAL1"]+np.arange(hdr["NAXIS1"])*hdr.get("CDELT1",hdr.get("CD1_1"))
    return wav,flux,interp1d(wav,flux,kind="linear",bounds_error=False,fill_value=np.nan)


def parse_orders(hdr):
    wat2=""
    i=1
    while f"WAT2_{i:03d}" in hdr: wat2+=f"{hdr[f'WAT2_{i:03d}']:<68}"; i+=1
    specs=re.findall(r'spec\d+\s*=\s*"([^"]+)"',wat2)
    return [{"ap":int(float(s.split()[0])),"w1":float(s.split()[3]),
             "dw":float(s.split()[4]),"nw":int(float(s.split()[5]))} for s in specs]


def wav_arr(o): return o["w1"]+np.arange(o["nw"])*o["dw"]


def _hdr_float(hdr, *keys):
    for k in keys:
        v = hdr.get(k)
        if v is None or isinstance(v, bool):
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if np.isfinite(f):
            return f
    return None


def _clock_str(v):
    """'HH:MM[:SS[.s]]' string, or seconds after midnight → 'HH:MM:SS.sss'."""
    if v is None or isinstance(v, bool):
        return None
    if isinstance(v, (int, float, np.integer, np.floating)):
        s = float(v)
        if not (0.0 <= s < 86400.0):
            return None
        h = int(s // 3600); m = int((s % 3600) // 60)
        return f"{h:02d}:{m:02d}:{s - 3600 * h - 60 * m:06.3f}"
    s = str(v).strip()
    return s if re.match(r"^\d{1,2}:\d{2}(:\d{2}(\.\d*)?)?$", s) else None


def obs_time_from_header(hdr):
    """Exposure start and mid-exposure UTC times.

    Tries, in order: DATE-OBS as ISO date-time; DATE-OBS date + TIME-OBS /
    UT / UT-START; UT-DATE + UT-START / UT-TIME (CarPy/MIKE); MJD-OBS / MJD;
    JD.  Each is taken as the exposure START; mid = start + EXPTIME/2.
    Returns (t_start, t_mid, source) or (None, None, reason).
    """
    t0, src = None, None
    try:
        dob = hdr.get("DATE-OBS")
        dob = str(dob).strip() if dob else ""
        if "T" in dob:
            t0, src = Time(dob, format="isot", scale="utc"), "DATE-OBS"
        for date_key, clock_keys in (
                ("DATE-OBS", ("TIME-OBS", "UT", "UT-START", "UTSTART", "UT-TIME")),
                ("UT-DATE",  ("UT-START", "UTSTART", "UT-TIME", "UT"))):
            if t0 is not None:
                break
            d = hdr.get(date_key)
            d = str(d).strip()[:10] if d else ""
            if not re.match(r"^\d{4}-\d{2}-\d{2}$", d):
                continue
            for ck in clock_keys:
                clk = _clock_str(hdr.get(ck))
                if clk:
                    t0, src = Time(f"{d}T{clk}", format="isot", scale="utc"), f"{date_key}+{ck}"
                    break
        if t0 is None:
            mjd = _hdr_float(hdr, "MJD-OBS", "MJD")
            if mjd is not None:
                t0, src = Time(mjd, format="mjd", scale="utc"), "MJD-OBS"
        if t0 is None:
            jd = _hdr_float(hdr, "JD")
            if jd is not None:
                t0, src = Time(jd, format="jd", scale="utc"), "JD"
    except Exception as e:
        return None, None, f"could not parse the observation time ({e})"
    if t0 is None:
        return None, None, ("no observation time in header (looked for DATE-OBS, "
                            "UT-DATE+UT-START/UT-TIME, MJD-OBS, JD)")
    exp = _hdr_float(hdr, "EXPTIME", "EXPOSURE") or 0.0
    return t0, t0 + 0.5 * exp * u.second, src + (" + EXPTIME/2" if exp else " (start)")


def obs_site_from_header(hdr):
    """EarthLocation from SITELAT/SITELONG (east-positive)/SITEALT, else
    OBSERVATORY (Las Campanas).  The site only enters through Earth's
    rotation, < 0.5 km/s."""
    try:
        lat, lon = hdr.get("SITELAT"), hdr.get("SITELONG")
        if lat is not None and lon is not None:
            alt = _hdr_float(hdr, "SITEALT") or 0.0
            return (EarthLocation(lat=Angle(lat, unit=u.deg), lon=Angle(lon, unit=u.deg),
                                  height=alt * u.m),
                    f"SITELAT/SITELONG ({hdr.get('SITENAME', 'site')})")
    except Exception:
        pass
    return OBSERVATORY, "default site (Las Campanas)"


def obs_coord_from_header(hdr):
    """Pointing as SkyCoord: RA-D/DEC-D (degrees), else RA/DEC (sexagesimal
    hours/degrees, or degrees if numeric).  EQUINOX other than 2000 is
    precessed to ICRS."""
    ra, dec = _hdr_float(hdr, "RA-D", "RA_DEG"), _hdr_float(hdr, "DEC-D", "DEC_DEG")
    if ra is not None and dec is not None:
        c, src = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs"), "RA-D/DEC-D"
    else:
        ra, dec = hdr.get("RA"), hdr.get("DEC")
        if ra is None or dec is None:
            return None, "no RA/DEC in header"
        if isinstance(ra, (int, float)) and isinstance(dec, (int, float)):
            c = SkyCoord(float(ra) * u.deg, float(dec) * u.deg, frame="icrs")
        else:
            c = SkyCoord(str(ra).strip(), str(dec).strip(),
                         unit=(u.hourangle, u.deg), frame="icrs")
        src = "RA/DEC"
    eq = _hdr_float(hdr, "EQUINOX", "EPOCH")
    if eq is not None and abs(eq - 2000.0) > 0.5:
        c = SkyCoord(c.ra, c.dec, frame=FK5(equinox=Time(eq, format="jyear"))).icrs
        src += f" (EQUINOX {eq:g} → ICRS)"
    return c, src


def _obs_meta_cached(fp_str, mtime):
    try:
        hdr = fits.getheader(fp_str)
    except Exception as e:
        return {"date": "—", "mjd": "—", "v_bary": float("nan"),
                "note": f"could not read header ({e})"}
    t0, tmid, tsrc = obs_time_from_header(hdr)
    out = {"date": t0.isot[:19] if t0 is not None else "—",
           "mjd": round(float(tmid.mjd), 5) if tmid is not None else "—",
           "v_bary": float("nan"), "note": tsrc}
    if tmid is None:
        return out
    coord, csrc = obs_coord_from_header(hdr)
    if coord is None:
        out["note"] = csrc
        return out
    site, ssrc = obs_site_from_header(hdr)
    try:
        out["v_bary"] = float(coord.radial_velocity_correction(
            kind="barycentric", obstime=tmid, location=site).to(u.km / u.s).value)
        out["note"] = f"time: {tsrc}; pointing: {csrc}; {ssrc}"
    except Exception as e:
        out["note"] = f"barycentric correction failed ({e})"
    return out


def obs_meta(fp):
    """{'date': start ISO UTC or '—', 'mjd': mid-exposure MJD or '—',
    'v_bary': km/s or NaN, 'note': where each came from / why it failed}."""
    p = Path(fp)
    try:
        mt = p.stat().st_mtime
    except OSError:
        mt = 0.0
    return dict(_obs_meta_cached(str(p), mt))


def fit_continuum(wav, flux, ivar=None, func="spline", order=3, maxiter=5,
                  low_sig=2.0, high_sig=1.0, knot_spacing=20, exclude=None,
                  scale=0.0, edge_pixels=0):
    """
    Fit continuum to a 1-D spectrum, mirroring alexmods `fit_continuum_lsq`.

    The spline path uses scipy.interpolate.LSQUnivariateSpline with explicit
    interior knots and iterative noise-weighted sigma clipping — the *exact*
    estimator alexmods/SMHR use, so HASHEM's normalized spectra are directly
    comparable to published SMHR-based results.

    alexmods reference (continuum.py :: fit_continuum_lsq):
        sig = (cont - y) * sqrt(ivar);  sig /= nanstd(sig)
        mask[sig > sigma_hi] = False     # below continuum  → absorption
        mask[sig < -sigma_lo] = False    # above continuum  → cosmics
    Our (high_sig, low_sig) map to alexmods (sigma_hi, sigma_lo) identically.

    Edge handling
    ─────────────
    scale       : exclude pixels below scale × 95th-pct flux (blaze rolloff).
    edge_pixels : hard-trim N pixels from each end before fitting.
    snip_edges  : continuum set to NaN beyond the first/last fitted pixel
                  (alexmods default) so blaze wings never get extrapolated.
    """
    from scipy import interpolate
    import numpy.polynomial.legendre as _leg
    import numpy.polynomial.chebyshev as _cheb

    x  = np.asarray(wav,  dtype=float)
    y  = np.asarray(flux, dtype=float)
    if ivar is None:
        ivar = np.ones_like(y)
    w  = np.asarray(ivar, dtype=float)
    w[~np.isfinite(w)] = 0.0

    n = len(x)
    _fallback = (np.nanmedian(y[np.isfinite(y)]) if np.isfinite(y).any() else 1.0)

    # ── Persistent good-pixel mask (alexmods uses a boolean mask) ─────────────
    mask = np.isfinite(x) & np.isfinite(y) & (y > 0) & (w >= 0)
    mask[np.abs(y) < 1e-6] = False            # alexmods: kill ~zero fluxes

    if edge_pixels > 0 and n > 2 * edge_pixels:
        mask[:edge_pixels]  = False
        mask[-edge_pixels:] = False

    if scale > 0.0 and mask.any():
        _peak = np.nanpercentile(y[mask], 95)
        if _peak > 0:
            mask &= ~(y < scale * _peak)

    if exclude:
        for _xmin, _xmax in exclude:
            mask[(x >= _xmin) & (x <= _xmax)] = False

    if mask.sum() < max(order + 2, 6):
        return np.full_like(x, _fallback)

    # Record edges for snip_edges (first/last pixel that survived masking)
    left  = int(np.where(mask)[0][0])
    right = int(np.where(mask)[0][-1])

    # ── Interior knots (alexmods initialize_knots) ────────────────────────────
    d0, d1   = x[mask].min(), x[mask].max()
    waverange = d1 - d0 - 2 * knot_spacing
    use_spline = func == "spline" and waverange > knot_spacing
    if use_spline:
        n_knots  = int(waverange // knot_spacing)
        min_knot = d0 + (waverange - n_knots * knot_spacing) / 2.0
        knots    = list(np.arange(min_knot, d1, knot_spacing))
    else:
        knots = []

    k_deg = int(np.clip(order, 1, 5))
    continuum = np.full_like(x, np.nanmedian(y[mask]))
    fcont = None

    for _ in range(maxiter):
        if mask.sum() < max(k_deg + 1, 6):
            break
        xm, ym, wm = x[mask], y[mask], w[mask]

        try:
            if use_spline:
                # Trim knots to stay strictly inside the current data range
                wmin, wmax = xm.min(), xm.max()
                kn = list(knots)
                while kn and kn[-1] >= wmax:
                    kn = kn[:-1]
                while kn and kn[0] <= wmin:
                    kn = kn[1:]
                if not kn:
                    use_spline = False
                    continue
                fcont = interpolate.LSQUnivariateSpline(
                    xm, ym, kn, w=wm, k=k_deg)
                continuum = fcont(x)
            elif func in ("poly", "polynomial") or not use_spline:
                c = np.polyfit(xm, ym, min(order, 5),
                               w=np.sqrt(np.clip(wm, 0, None)))
                continuum = np.polyval(c, x)
            elif func in ("leg", "legendre"):
                c = _leg.legfit(xm, ym, order, w=np.sqrt(np.clip(wm, 0, None)))
                continuum = _leg.legval(x, c)
            else:  # chebyshev
                c = _cheb.chebfit(xm, ym, order, w=np.sqrt(np.clip(wm, 0, None)))
                continuum = _cheb.chebval(x, c)
        except Exception:
            break

        # ── Noise-weighted sigma clipping (alexmods convention) ───────────────
        # sig = (cont - y) * sqrt(ivar);   sig /= nanstd(sig)
        sig = (continuum - y) * np.sqrt(np.clip(w, 0, None))
        finite = np.isfinite(sig)
        if finite.sum() < 5:
            break
        _std = np.nanstd(sig[finite])
        if _std == 0:
            break
        sig = sig / _std
        new_mask = mask.copy()
        new_mask[sig > high_sig] = False      # flux below continuum (absorption)
        new_mask[sig < -low_sig] = False      # flux above continuum (cosmics)
        if np.array_equal(new_mask, mask):
            break
        mask = new_mask

    # Final evaluation
    if use_spline and fcont is not None:
        continuum = fcont(x)
    continuum = np.where(np.isfinite(continuum) & (continuum > 0),
                         continuum, _fallback)

    # ── snip_edges (alexmods default) ─────────────────────────────────────────
    continuum[:left]      = np.nan
    continuum[right + 1:] = np.nan
    return continuum


def normed_flux(o, data, i, func="spline", order=2, maxiter=5,
                low_sig=2.0, high_sig=1.0, knot_spacing=20):
    """Return the continuum-normalized OBJECT spectrum (band 2) for one order,
    exactly as the SMHR reference does: flux=band2, ivar=band3^-2, fit spline,
    divide. Falls back to band 6 (object/normed-flat) if the fit fails."""
    w   = wav_arr(o)
    obj = data[1, i, :]            # band 2 = object spectrum
    nz  = data[2, i, :]            # band 3 = noise spectrum
    ivar = np.where(nz > 0, nz**-2.0, 0.0)
    try:
        cont = fit_continuum(w, obj, ivar=ivar, func=func, order=order,
                             maxiter=maxiter, low_sig=low_sig,
                             high_sig=high_sig, knot_spacing=knot_spacing)
        fn = obj / np.where(cont > 0, cont, np.nan)
        if np.isfinite(fn).sum() < 10:
            return data[6, i, :]
        return fn
    except Exception:
        return data[6, i, :]


def joint_ccf(ap_order,ord_dict,data,rv_grid,anchor=None):
    _,_,interp_t=load_template(_TEMPLATE["path"])
    ccf=np.zeros(len(rv_grid)); tw=0.0
    rv_mid=rv_grid[len(rv_grid)//2] if anchor is None else anchor
    for i,ap in enumerate(ap_order):
        o=ord_dict[ap]; w=wav_arr(o); fl=normed_flux(o,data,i)
        war=w/(1+rv_mid/C_KMS)
        if not any(((war>=r0)&(war<=r1)).sum()>=10 for r0,r1,_ in CCF_REGIONS): continue
        oc=np.zeros(len(rv_grid))
        for k,rv in enumerate(rv_grid):
            wr=w/(1+rv/C_KMS); in_r=np.zeros(len(w),dtype=bool)
            for r0,r1,_ in CCF_REGIONS: in_r|=(wr>=r0)&(wr<=r1)
            if in_r.sum()<10: continue
            f,t=fl[in_r],interp_t(wr[in_r]); v=np.isfinite(f)&np.isfinite(t)&(f>0)&(t>0)
            if v.sum()<10: continue
            oc[k]=np.dot(f[v]/np.median(f[v])-1,t[v]-1)
        if anchor is not None:
            near=np.abs(rv_grid-anchor)<=30
            wt=float(np.clip(oc[near].max(),0,None)) if near.any() else 0.0
        else:
            pk=oc.max(); wt=(1.0 if pk>1e-10 else 0.0)
            if wt: oc/=pk
        if wt>0: ccf+=oc*(wt if anchor else 1.0); tw+=wt
    return ccf/tw if tw>0 else ccf


def para_peak(rg,ccf):
    pi=np.argmax(ccf)
    if 2<=pi<=len(ccf)-3:
        x,y=rg[pi-2:pi+3],ccf[pi-2:pi+3]; c=np.polyfit(x,y,2); return -c[1]/(2*c[0])
    return rg[pi]

def measure_rv(fp,ap_order,ord_dict,data):
    rc=np.arange(-500,500,2.); cc=joint_ccf(ap_order,ord_dict,data,rc)
    ra=para_peak(rc,cc); rf=np.arange(ra-50,ra+50+0.5,0.5)
    cf=joint_ccf(ap_order,ord_dict,data,rf,anchor=ra); rv=para_peak(rf,cf)
    return rv,float(cf.max()),rf,cf


def _tukey(n, alpha=0.2):
    """Tukey (tapered cosine) window; alpha in [0,1] = fraction tapered."""
    try:
        from scipy.signal.windows import tukey as _sp_tukey
        return _sp_tukey(n, alpha)
    except Exception:
        if alpha <= 0:
            return np.ones(n)
        w = np.ones(n)
        edge = int(alpha * (n - 1) / 2.0)
        if edge < 1:
            return w
        k = np.arange(edge)
        taper = 0.5 * (1 + np.cos(np.pi * (k / edge - 1)))
        w[:edge] = taper
        w[-edge:] = taper[::-1]
        return w


def _antisym_centroid(vel, ccf, pidx, win_k, search_pix=1.5, nsub=61):
    """
    Tonry & Davis (1979) antisymmetric-component peak centroid.

    The true velocity is the center about which the CCF is maximally SYMMETRIC:
    a perfect, noise-free correlation peak is an even function of the shift, so
    the odd (antisymmetric) part about the true peak is pure noise.  We search a
    sub-pixel offset δ around the integer peak and pick the δ that minimizes the
    RMS of the antisymmetric part
        a_k(δ) = ½ [ c(v_pk+δ+kΔ) − c(v_pk+δ−kΔ) ],  k = 1 … win_k
    over the whole peak window.  Because it uses the entire peak shape rather
    than the three pixels a parabola sees, it is markedly more robust to noise
    at the very top of the CCF (which is exactly where a 3-point parabola is
    weakest).  Returns (rv, peak_height, curvature); curvature>0 at a real max
    and feeds the FWHM used for the velocity error.
    """
    dv     = vel[1] - vel[0]
    v_pk   = vel[pidx]
    ks     = (np.arange(1, win_k + 1) * dv)[None, :]                 # (1, win_k)
    deltas = (np.linspace(-search_pix, search_pix, nsub) * dv)[:, None]  # (nsub,1)
    centers = v_pk + deltas
    c_plus  = np.interp((centers + ks).ravel(), vel, ccf).reshape(nsub, win_k)
    c_minus = np.interp((centers - ks).ravel(), vel, ccf).reshape(nsub, win_k)
    asym    = np.mean((0.5 * (c_plus - c_minus)) ** 2, axis=1)       # (nsub,)

    j = int(np.argmin(asym))
    if 1 <= j <= len(asym) - 2:                    # parabolic refine of the min
        y0, y1, y2 = asym[j - 1], asym[j], asym[j + 1]
        den  = y0 - 2 * y1 + y2
        frac = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        d_best = float(deltas[j, 0] + frac * (deltas[1, 0] - deltas[0, 0]))
    else:
        d_best = float(deltas[j, 0])
    rv = v_pk + d_best

    peak_h = float(np.interp(rv, vel, ccf))
    c_l = float(np.interp(rv - dv, vel, ccf))
    c_r = float(np.interp(rv + dv, vel, ccf))
    curv = -(c_l - 2.0 * peak_h + c_r)             # >0 for a maximum
    return rv, peak_h, curv

def tonry_davis_ccf(wav, flux, template_interp, vmin=-500.0, vmax=500.0,
                    apodize=0.2):
    """
    Tonry & Davis (1979) FFT cross-correlation on a log-lambda grid.

    Parameters
    ----------
    wav, flux : 1-D arrays
        Object spectrum (continuum-flattened band is fine; we re-flatten).
    template_interp : callable
        interp1d over the rest-frame template flux vs wavelength.
    vmin, vmax : float
        Velocity search window (km/s).
    apodize : float
        Tukey taper fraction applied to both spectra before the FFT.

    Returns
    -------
    (rv, rv_err, r_value, vel_axis, ccf)  or  None if the fit is impossible.
        rv, rv_err in km/s; r_value is the Tonry-Davis detection statistic.
    """
    wav  = np.asarray(wav, float)
    flux = np.asarray(flux, float)
    good = np.isfinite(wav) & np.isfinite(flux) & (flux > 0)
    if good.sum() < 20:
        return None
    wav, flux = wav[good], flux[good]

    # Sort and remove non-increasing (overlapping orders) wavelengths
    srt = np.argsort(wav)
    wav, flux = wav[srt], flux[srt]
    keep = np.concatenate([[True], np.diff(wav) > 0])
    wav, flux = wav[keep], flux[keep]
    if len(wav) < 20:
        return None

    # ── Uniform ln(λ) grid: velocity shift ↔ constant pixel lag ──────────────
    lnw   = np.log(wav)
    dln   = np.median(np.diff(lnw))
    if not np.isfinite(dln) or dln <= 0:
        return None
    n_pix = int((lnw[-1] - lnw[0]) / dln)
    if n_pix < 32:
        return None
    lng   = lnw[0] + np.arange(n_pix) * dln
    wavg  = np.exp(lng)

    # Interpolate object + template onto the grid
    f_obj  = np.interp(wavg, wav, flux, left=np.nan, right=np.nan)
    f_tmpl = template_interp(wavg)
    m = np.isfinite(f_obj) & np.isfinite(f_tmpl) & (f_tmpl > 0)
    if m.sum() < 32:
        return None

    # Continuum-flatten both to fluctuate about zero (divide by median, −1)
    f_obj  = f_obj  / np.nanmedian(f_obj[m])  - 1.0
    f_tmpl = f_tmpl / np.nanmedian(f_tmpl[m]) - 1.0
    f_obj[~m]  = 0.0
    f_tmpl[~m] = 0.0

    # Apodize the ends to suppress FFT ringing
    win = _tukey(n_pix, apodize)
    f_obj  *= win
    f_tmpl *= win

    # Normalize to unit RMS so the CCF peak ≈ correlation coefficient
    so  = np.sqrt(np.mean(f_obj  ** 2))
    stp = np.sqrt(np.mean(f_tmpl ** 2))
    if so <= 0 or stp <= 0:
        return None
    f_obj  /= so
    f_tmpl /= stp

    # ── FFT cross-correlation ────────────────────────────────────────────────
    F_o = np.fft.rfft(f_obj)
    F_t = np.fft.rfft(f_tmpl)
    ccf = np.fft.irfft(F_o * np.conj(F_t), n=n_pix)
    ccf = np.fft.fftshift(ccf) / n_pix

    lags = np.arange(n_pix) - n_pix // 2
    vel  = lags * dln * C_KMS

    sel = (vel >= vmin) & (vel <= vmax)
    if sel.sum() < 5:
        return None
    vsel = vel[sel]
    csel = ccf[sel]

    pi = int(np.argmax(csel))
    dv_step = vsel[1] - vsel[0]
    # peak index in the FULL ccf array (need room for the symmetry window)
    pidx = int(np.argmin(np.abs(vel - vsel[pi])))
    half = min(pidx, n_pix - 1 - pidx)
    win_k = min(half, max(8, int(60.0 / abs(dv_step)) if dv_step else 8))

    if win_k >= 3:
        # ── Tonry-Davis antisymmetric-component centroid ──
        rv, peak_h, curv = _antisym_centroid(vel, ccf, pidx, win_k)
        # Noise proxy for r: RMS of the antisymmetric part about the INTEGER
        # peak — kept fixed (not the minimized centroid) so r is not biased
        # high by "measuring what we just minimized".
        ks = np.arange(1, win_k + 1)
        a  = (ccf[pidx + ks] - ccf[pidx - ks]) / 2.0
        sigma_a = np.sqrt(np.mean(a ** 2))
    else:
        # Fallback: 3-point parabola when the window is too small
        if 1 <= pi <= len(csel) - 2:
            y0, y1, y2 = csel[pi - 1], csel[pi], csel[pi + 1]
            denom = (y0 - 2 * y1 + y2)
            frac  = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
            rv     = vsel[pi] + frac * dv_step
            peak_h = y1 - 0.25 * (y0 - y2) * frac
            curv   = -denom
        else:
            rv, peak_h, curv = vsel[pi], csel[pi], 0.0
        sigma_a = np.std(csel)

    h = peak_h
    r_value = float(h / (np.sqrt(2.0) * sigma_a)) if sigma_a > 0 else 0.0

    # Velocity error (Tonry-Davis): σ_v ≈ (3/8) · w / (1 + r)
    if curv > 0 and peak_h > 0:
        sigma_pix = np.sqrt(abs(peak_h / curv))
        fwhm = 2.3548 * sigma_pix * dv_step
    else:
        fwhm = 3.0 * abs(dv_step)
    rv_err = float((3.0 / 8.0) * abs(fwhm) / (1.0 + max(r_value, 0.0)))

    return rv, rv_err, r_value, vel, ccf

RV_SEARCH_VMAX = 500.0


def measure_rv_region_td(ap_order, ord_dict, data, region, template_interp,
                         vmax=RV_SEARCH_VMAX):
    """
    Run one Tonry-Davis CCF for a single CCF region by concatenating the
    pixels of every order that overlaps that region (rest frame).
    Returns (rv, rv_err, r_value, vel_axis, ccf) or None.
    """
    pieces = ccf_region_pixels(ap_order, ord_dict, data, region)
    if not pieces:
        return None
    return tonry_davis_ccf(np.concatenate([p[1] for p in pieces]),
                           np.concatenate([p[2] for p in pieces]),
                           template_interp, vmin=-vmax, vmax=vmax)


def ccf_region_pixels(ap_order, ord_dict, data, region):
    """The pixels one CCF region correlates: every order with >= 10 pixels
    within the region +/- 5 A (OBSERVED wavelengths), continuum-normalized
    with normed_flux() defaults.  Returns [(ap, wav, flux), ...].  Shared by
    measure_rv_region_td and the RV tab's spectrum check, so the plot shows
    exactly the CCF input (tonry_davis_ccf then sorts, resamples to ln(lambda)
    and divides by the median)."""
