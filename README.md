# RV-correction

Radial velocities of Magellan/MIKE echelle spectra by cross-correlation with a
rest-frame template, with a figure of each spectrum against the template.
It is the RV measurement of [HASHEM](https://github.com/Mohammad-Mardini/HASHEM)
as a single script, and gives the same RVs as the HASHEM app.

Author: Mohammad K. H. Mardini

## Method

Each CCF region (Ca II triplet 8450-8750 A, Halpha 6510-6610 A, Mg b
5100-5200 A, Hbeta 4810-4910 A, and 4290-4390 A) is continuum-normalized,
resampled to log wavelength, and cross-correlated with the template by FFT.
The RV comes from the Tonry & Davis (1979) antisymmetric-component centroid,
with their r statistic and error, searched within +/-500 km/s. The adopted RV
is that of the chosen region (default Ca II triplet), or of the highest-r
region if that one is not detected. When two or more regions have r >= 4, the
quoted error is the larger of the formal error and their scatter.
Barycentric corrections use the observation time, pointing and site in each
header.

## Installation

```
git clone https://github.com/Mohammad-Mardini/RV-correction.git
cd RV-correction
pip install -r requirements.txt
```

## Usage

Edit the settings at the top of `measure_rv.py` once:

```python
SPECTRA  = "*_multi.fits"       # spectra to measure (a glob pattern or a list)
TEMPLATE = "hd122563.fits"      # rest-frame RV template (1-D FITS)
REGION   = "Ca II IR"           # "Ca II IR", "Halpha", "Mg b", "Hbeta", "Ca H&K" or "auto"
OUT_CSV  = "rv_results.csv"     # results table
PLOT_DIR = "rv_plots"           # one figure per spectrum
PLOTS    = True                 # False: measure only
```

then run

```
python measure_rv.py
```

Command-line arguments override the settings:

```
python measure_rv.py --template path/to/template.fits --region auto path/to/*_multi.fits
```

## Input

* MIKE `*_multi.fits` spectra from CarPy (object, noise and S/N bands; linear
  WAT2 wavelength solutions). Co-added spectra (`*_coadd_*`) are skipped.
* A rest-frame 1-D template in FITS with a linear CRVAL1/CDELT1 solution.

## Output

* A table on screen and `rv_results.csv`: RV and error (topocentric), adopted
  region, Tonry-Davis r, v_bary and v_helio, the RV and r of every region,
  the observation date, and the template name. v_bary is left empty when the
  header lacks what it needs.
* `rv_plots/<spectrum>_rv.png`: for each region, the CCF (top) and the
  spectrum the CCF used with the template at the adopted RV (bottom).

## Citation

Please cite using `CITATION.cff` (GitHub shows a "Cite this repository" button).

## Licence

MIT; see `LICENSE`.
