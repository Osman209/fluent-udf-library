"""
wave_theory_check.py

Independent numpy re-implementation of the Airy and Stokes second order
formulas (Dean & Dalrymple) used in common/wave_theory.h. Reads the
reference samples written by tests/test_wave_theory.c and compares.
Any disagreement above 1e-9 fails.

Run:  python3 validation/wave_theory_check.py tests/out
"""
import sys, os, csv, math
import numpy as np
from scipy.optimize import brentq

out = sys.argv[1] if len(sys.argv) > 1 else "tests/out"
path = os.path.join(out, "wave_ref.csv")

with open(path) as f:
    header = f.readline().strip("# \n")
    params = dict(kv.split("=") for kv in header.split())
    rows = list(csv.DictReader(f))

H, T, d, g = (float(params[k]) for k in ("H", "T", "d", "g"))
A = H / 2
w = 2 * math.pi / T

# dispersion by an independent method (bracketing root finder)
k = brentq(lambda kk: g * kk * math.tanh(kk * d) - w * w, 1e-6, 1e3)
assert abs(k - float(params["k"])) / k < 1e-9, "wavenumber mismatch with C"

def airy(x, y, t):
    th = k * x - w * t
    eta = A * np.cos(th)
    yy = min(y, eta)
    u = A * w * np.cosh(k * (yy + d)) / np.sinh(k * d) * np.cos(th)
    v = A * w * np.sinh(k * (yy + d)) / np.sinh(k * d) * np.sin(th)
    return eta, u, v

def stokes2(x, y, t):
    th = k * x - w * t
    sh, ch = np.sinh(k * d), np.cosh(k * d)
    eta = A * np.cos(th) + (k * A**2 / 4) * ch / sh**3 * (2 + np.cosh(2 * k * d)) * np.cos(2 * th)
    yy = min(y, eta)
    u = (A * w * np.cosh(k * (yy + d)) / sh * np.cos(th)
         + 0.75 * A**2 * w * k * np.cosh(2 * k * (yy + d)) / sh**4 * np.cos(2 * th))
    v = (A * w * np.sinh(k * (yy + d)) / sh * np.sin(th)
         + 0.75 * A**2 * w * k * np.sinh(2 * k * (yy + d)) / sh**4 * np.sin(2 * th))
    return eta, u, v

worst = 0.0
for r in rows:
    fn = airy if r["theory"] == "airy" else stokes2
    eta, u, v = fn(float(r["x"]), float(r["y"]), float(r["t"]))
    err = max(abs(eta - float(r["eta"])), abs(u - float(r["u"])), abs(v - float(r["v"])))
    worst = max(worst, err)

print(f"wave_theory_check: {len(rows)} samples, worst |C - numpy| = {worst:.2e}")
if worst > 1e-9:
    print("FAIL")
    sys.exit(1)
print("ok: C wave theory matches independent numpy implementation")
