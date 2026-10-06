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

