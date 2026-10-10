# RV-correction

![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)

**Radial velocities for multi-epoch Magellan/MIKE echelle spectra, and rest-frame
copies of every spectrum ready for co-addition.**

One script, `measure_rv.py`. Point it at your `*_multi.fits` files and an RV
template, run it, and get:

- a table of radial velocities with uncertainties and barycentric corrections,
- a diagnostic figure per spectrum (CCF and spectrum versus template), and
- each spectrum shifted to the rest frame by its own RV, as FITS.

It is the RV step of [HASHEM](https://github.com/Mohammad-Mardini/HASHEM), packaged
to run on its own. It gives the same velocities as the HASHEM app.

---

## Quick start

```bash
git clone https://github.com/Mohammad-Mardini/RV-correction.git
cd RV-correction
pip install -r requirements.txt
```

Open `measure_rv.py`, set the template path (and anything else) in the settings
block at the top, then run it from the folder that holds your spectra:

```bash
python measure_rv.py
```

That is all. Results appear on screen and in the files listed under
[Output](#output).

---

## Settings

All settings sit at the top of `measure_rv.py`. Relative paths are taken from
the folder you run the script in.

| Setting | Default | What it does |
|---|---|---|
| `SPECTRA` | `"*_multi.fits"` | Spectra to measure: a glob pattern or a list of patterns |
| `TEMPLATE` | `"hd122563.fits"` | Rest-frame RV template (1-D FITS, linear `CRVAL1`/`CDELT1`) |
| `REGION` | `"Ca II IR"` | Region whose RV is adopted: `"Ca II IR"`, `"Halpha"`, `"Mg b"`, `"Hbeta"`, `"Ca H&K"`, or `"auto"` |
| `OUT_CSV` | `"rv_results.csv"` | Results table |
| `PLOT_DIR` | `"rv_plots"` | Folder for the diagnostic figures |
| `PLOTS` | `True` | `False` to skip the figures |
| `SAVE_SHIFTED` | `True` | Write rest-frame copies of the spectra |
| `SHIFTED_DIR` | `"rv_shifted"` | Folder for the rest-frame copies |
| `OBSERVATORY` | `"Magellan (Las Campanas; MIKE)"` | Site for the barycentric correction (see [Observatory](#observatory)) |
| `USE_HEADER_SITE` | `True` | Use each header's own site keywords when present |

Every setting can also be given on the command line, which overrides the file:

```bash
python measure_rv.py --template ref/hd122563.fits --region auto data/*_multi.fits
python measure_rv.py --help
```

---

## Output

| File | Contents |
|---|---|
| `rv_results.csv` | One row per spectrum: RV and error (topocentric, km/s), adopted region, Tonry-Davis *r*, `v_bary`, `v_helio`, the RV and *r* of every region, observation date and MJD, template, site, and the name of the rest-frame copy |
| `rv_plots/<name>_rv.png` | Top: the CCF of each detected region, with the adopted region marked and a warning if the peak sits at the edge of the search window. Bottom: the spectrum the CCF actually used, with the template drawn at the adopted RV |
| `rv_shifted/<name>__RV_shifted_<arm>_multi.fits` | The spectrum in the rest frame. Only the wavelength solution changes; the pixel data are copied unchanged |

Headers of the rest-frame copies record what was done:

| Keyword | Meaning |
|---|---|
| `RV_KMS` | RV the spectrum was shifted by (km/s) |
| `RV_ERR` | Its uncertainty (km/s) |
| `RV_REGION` | Region the RV came from |
| `RV_TMPL` | Template the RV was measured with |
| `RV_MEAS` | RV as measured |
| `V_BARY`, `V_HELIO` | Barycentric correction and corrected RV (km/s); omitted when the header lacks what the correction needs |

Co-added spectra (`*_coadd_*`) and files that are already rest-frame copies
(`*_RV_shifted_*`) are skipped as input, so running the script twice in the
same folder does not shift anything twice.

---

## How the RV is measured

For each spectrum and each CCF region:

1. **Select pixels.** Every echelle order with at least 10 pixels within the
   region (plus 5 Å on each side) contributes.
2. **Normalize.** Each order is continuum-normalized with a spline (20 Å knot
   spacing) fitted with iterative, noise-weighted sigma clipping.
3. **Cross-correlate.** Spectrum and template are resampled to a common uniform
   ln λ grid, divided by their median, tapered with a cosine (Tukey) window, and
   cross-correlated by FFT. A constant lag on this grid is a constant velocity.
4. **Locate the peak.** Within ±500 km/s, the velocity comes from the
   antisymmetric-component centroid of the CCF peak (Tonry & Davis 1979), with
   their statistic *r* = *h* / (√2 σ<sub>a</sub>) and error σ<sub>v</sub> = (3/8) *w* / (1 + *r*).

The velocities are combined as follows:

- **Adopted RV:** that of the region set by `REGION`. If that region is not
  detected, or `REGION = "auto"`, the region with the highest *r* is used. If no
  region is detected, a joint CCF over all regions is used.
- **Uncertainty:** when two or more regions have *r* ≥ 4, the quoted error is the
  larger of the formal error and the scatter of their velocities, so disagreement
  between regions is not hidden.
- **Barycentric correction:** computed with astropy at mid-exposure, from the
  observation time (`DATE-OBS`, or MIKE's `UT-DATE` with `UT-START`/`UT-TIME`,
  or `MJD-OBS`), the pointing (`RA-D`/`DEC-D` or `RA`/`DEC`), and the site.
  `v_helio` = RV + `v_bary`.
- **Rest frame:** λ<sub>rest</sub> = λ<sub>obs</sub> / (1 + RV/*c*), applied to
  each order's wavelength solution.

The RVs in the table and in `RV_KMS` are topocentric: they are what the
spectrum must be shifted by, because each spectrum's wavelengths are in the
frame of its own night.

