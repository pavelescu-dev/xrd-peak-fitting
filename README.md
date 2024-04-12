# xrd-peak-fitting

Small helper script I use for fitting XRD peaks with pseudo-Voigt profiles.
Written for my own convenience — not production quality.

## Install

```
pip install -r requirements.txt
```

## Usage

```
python fit.py --input sample.xy --peaks 2theta=28.4,32.1,45.7 --out fit.json
```

Outputs peak positions, FWHM, integrated intensity, and a residuals plot.

## Why

The in-house software we have at the lab is licensed per-seat and the CLI is
clunky. I just wanted something I can run on my laptop while on the train.

## References

- Young, R.A. (1993) *The Rietveld Method*
- scipy.optimize.curve_fit docs
