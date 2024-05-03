#!/usr/bin/env python3
"""
XRD peak fitting with pseudo-Voigt profiles.
Inputs a 2-column file (2theta, intensity) and a list of peak centers.
"""
import argparse
import json
import numpy as np
from scipy.optimize import curve_fit


def pseudo_voigt(x, amp, center, fwhm, eta):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gamma = fwhm / 2.0
    gaussian = np.exp(-((x - center) ** 2) / (2.0 * sigma ** 2))
    lorentzian = 1.0 / (1.0 + ((x - center) / gamma) ** 2)
    return amp * (eta * lorentzian + (1.0 - eta) * gaussian)


def multi_peak(x, *params):
    n = len(params) // 4
    y = np.zeros_like(x)
    for i in range(n):
        amp, center, fwhm, eta = params[4 * i: 4 * i + 4]
        y = y + pseudo_voigt(x, amp, center, fwhm, eta)
    return y


def load_xy(path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def fit(x, y, centers, fwhm_init=0.2):
    p0 = []
    for c in centers:
        # amp ~ local maximum near c
        mask = np.abs(x - c) < fwhm_init * 3
        amp_guess = y[mask].max() if mask.any() else y.max()
        p0 += [amp_guess, c, fwhm_init, 0.5]
    popt, pcov = curve_fit(multi_peak, x, y, p0=p0, maxfev=20000)
    result = []
    for i in range(len(centers)):
        amp, ctr, fwhm, eta = popt[4 * i: 4 * i + 4]
        # integrated intensity approx (trapezoid of single component)
        yc = pseudo_voigt(x, amp, ctr, fwhm, eta)
        integ = float(np.trapz(yc, x))
        result.append({
            "center_2theta": float(ctr),
            "fwhm": float(fwhm),
            "amp": float(amp),
            "eta": float(eta),
            "integrated_intensity": integ,
        })
    residuals = y - multi_peak(x, *popt)
    return result, residuals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--peaks", required=True,
                    help="comma-separated 2theta centers, e.g. 28.4,32.1")
    ap.add_argument("--out", default="fit.json")
    args = ap.parse_args()

    x, y = load_xy(args.input)
    centers = [float(v) for v in args.peaks.replace("2theta=", "").split(",")]
    peaks, resid = fit(x, y, centers)
    out = {
        "input": args.input,
        "centers": centers,
        "peaks": peaks,
        "rms_residual": float(np.sqrt(np.mean(resid ** 2))),
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {args.out} -- rms={out['rms_residual']:.4g}")


if __name__ == "__main__":
    main()
