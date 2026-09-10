"""
motion_reference.py

Reference solutions for the prescribed-motion UDFs, so a Fluent run can
be compared against something exact.

1. cg_motion_translate.c: half-cosine velocity ramp. Checks that the
   displacement at the end of the ramp is V * T_RAMP / 2 (integral of
   the ramp).

2. cg_motion_oscillate.c: shows why the ramp is applied to the
   amplitude and the exact derivative returned. Integrates both
   variants numerically and reports the mean-position drift of the
   velocity-ramp variant (the one that was NOT used) versus zero drift
   for the amplitude-ramp variant (the one in the UDF).

Run:  python3 validation/motion_reference.py
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid

V, TR = 10.0, 1.0
t = np.linspace(0, 3, 30001)
r = np.where(t < TR, 0.5 * (1 - np.cos(np.pi * t / TR)), 1.0)
x = cumulative_trapezoid(V * r, t, initial=0)
x_end = np.interp(TR, t, x)
print(f"translate: x(T_RAMP) = {x_end:.6f}, expected V*T_RAMP/2 = {V*TR/2:.6f}")
assert abs(x_end - V * TR / 2) < 1e-4

H0, f, TR2 = 0.05, 1.0, 0.5
w = 2 * np.pi * f
r2 = np.where(t < TR2, 0.5 * (1 - np.cos(np.pi * t / TR2)), 1.0)
dr2 = np.where(t < TR2, 0.5 * (np.pi / TR2) * np.sin(np.pi * t / TR2), 0.0)
# variant A (used): exact derivative of r(t) H0 sin(wt)
vA = dr2 * H0 * np.sin(w * t) + r2 * H0 * w * np.cos(w * t)
hA = cumulative_trapezoid(vA, t, initial=0)
# variant B (not used): ramp multiplies the velocity only
vB = r2 * H0 * w * np.cos(w * t)
hB = cumulative_trapezoid(vB, t, initial=0)
mask = t > 2.0
driftA = hA[mask].mean(); driftB = hB[mask].mean()
errA = np.max(np.abs(hA[mask] - H0 * np.sin(w * t[mask])))
print(f"oscillate: amplitude-ramp variant follows H0 sin(wt) after ramp, max error {errA:.2e} m, mean drift {driftA:.2e} m")
print(f"oscillate: velocity-ramp variant (not used) has mean drift {driftB:.4f} m = {100*driftB/H0:.1f} % of H0")
assert errA < 1e-6
print("ok: motion reference checks passed")
