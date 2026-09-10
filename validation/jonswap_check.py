"""
jonswap_check.py

Rebuilds the JONSWAP component list (amplitude, wavenumber, omega,
phase) with numpy using the same LCG seed, and compares with the list
written by tests/test_wave_theory.c (jonswap_ref.csv). Also checks the
irregular-wave sample in irregular_ref.csv, and reports Hs recovered as
4 sqrt(m0).

The same script can read the wave_components.csv that the Fluent UDF
writes, to confirm the sea state a run actually used:
    python3 validation/jonswap_check.py <dir containing the csv> --fluent

Run:  python3 validation/jonswap_check.py tests/out
"""
import sys, os, csv, math
import numpy as np
from scipy.optimize import brentq

out = sys.argv[1] if len(sys.argv) > 1 else "tests/out"
fluent_mode = "--fluent" in sys.argv
fname = "wave_components.csv" if fluent_mode else "jonswap_ref.csv"
path = os.path.join(out, fname)

with open(path) as f:
    header = f.readline().strip("# \n")
    params = dict(kv.split("=") for kv in header.split())
    rows = list(csv.DictReader(f))

Hs, Tp, gam = float(params["Hs"]), float(params["Tp"]), float(params["gamma"])
N = int(params["N"]); fmin, fmax = float(params["fmin"]), float(params["fmax"])
d = float(params["d"]); seed = int(params["seed"]); g = 9.81
fp = 1 / Tp

def jonswap(f):
    sigma = np.where(f <= fp, 0.07, 0.09)
    r = np.exp(-(f / fp - 1) ** 2 / (2 * sigma ** 2))
    betaJ = 0.0624 * (1.094 - 0.01915 * math.log(gam)) / (0.230 + 0.0336 * gam - 0.185 / (1.9 + gam))
    return betaJ * Hs ** 2 * fp ** 4 * f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gam ** r

# same LCG as wave_theory.h
state = seed
def lcg():
    global state
    state = (1664525 * state + 1013904223) % 2 ** 32
    return state / 2 ** 32

df = (fmax - fmin) / N
f = fmin + (np.arange(N) + 0.5) * df
S = jonswap(f)
a = np.sqrt(2 * S * df)
w = 2 * np.pi * f
k = np.array([brentq(lambda kk: g * kk * math.tanh(kk * d) - ww * ww, 1e-6, 1e4) for ww in w])
p = np.array([2 * math.pi * lcg() for _ in range(N)])
m0 = np.sum(S * df)

ca = np.array([float(r["amplitude"]) for r in rows])
ck = np.array([float(r["wavenumber"]) for r in rows])
cw = np.array([float(r["omega"]) for r in rows])
cp = np.array([float(r["phase"]) for r in rows])

tol = 1e-7 if fluent_mode else 1e-9   # the UDF writes with %.10g
err = max(np.max(np.abs(a - ca)), np.max(np.abs(k - ck)), np.max(np.abs(w - cw)), np.max(np.abs(p - cp)))
print(f"jonswap_check: N={N}, worst component mismatch = {err:.2e}, 4 sqrt(m0) = {4*math.sqrt(m0):.4f} (Hs={Hs})")
if err > tol:
    print("FAIL: component list differs from C"); sys.exit(1)

if not fluent_mode:
    with open(os.path.join(out, "irregular_ref.csv")) as fh:
        r = list(csv.DictReader(fh))[0]
    x, y, t = float(r["x"]), float(r["y"]), float(r["t"])
    eta = np.sum(a * np.cos(k * x - w * t + p))
    yy = min(y, eta)
    u = np.sum(a * w * np.cosh(k * (yy + d)) / np.sinh(k * d) * np.cos(k * x - w * t + p))
    v = np.sum(a * w * np.sinh(k * (yy + d)) / np.sinh(k * d) * np.sin(k * x - w * t + p))
    e2 = max(abs(eta - float(r["eta"])), abs(u - float(r["u"])), abs(v - float(r["v"])))
    print(f"irregular sample mismatch = {e2:.2e}")
    if e2 > 1e-9:
        print("FAIL"); sys.exit(1)
print("ok: JONSWAP components and irregular kinematics match")
