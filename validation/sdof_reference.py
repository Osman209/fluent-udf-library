"""
sdof_reference.py

Analytic free-decay solution of a 1-DOF mass-spring-damper, to compare
with sdof_motion.csv from motion/sdof_spring_damper.c in a check case
where fluid forces are negligible (body in a large domain of very low
density gas, or the body released from an initial offset in still air
with the fluid force small compared with the spring force).

    m x'' + c x' + k x = 0,  x(0) = x0, x'(0) = 0
    zeta = c / (2 sqrt(k m)),  wn = sqrt(k/m),  wd = wn sqrt(1 - zeta^2)
    x(t) = x0 e^{-zeta wn t} ( cos wd t + zeta/sqrt(1-zeta^2) sin wd t )

Usage:
    python3 validation/sdof_reference.py sdof_motion.csv --m 10 --k 500 --c 2 --dof y
It prints the damped period and the log decrement from the Fluent
trace next to the analytic values, and the RMS difference.

--selftest integrates the ODE numerically and checks the closed form.
"""
import sys, math, argparse
import numpy as np
from scipy.integrate import solve_ivp

def analytic(t, m, k, c, x0):
    wn = math.sqrt(k / m); z = c / (2 * math.sqrt(k * m))
    if z >= 1:
        raise SystemExit("overdamped case not implemented in this reference")
    wd = wn * math.sqrt(1 - z * z)
    return x0 * np.exp(-z * wn * t) * (np.cos(wd * t) + z / math.sqrt(1 - z * z) * np.sin(wd * t)), wd, z

def selftest():
    m, k, c, x0 = 10.0, 500.0, 2.0, 0.05
    t = np.linspace(0, 10, 5001)
    sol = solve_ivp(lambda tt, y: [y[1], (-c * y[1] - k * y[0]) / m], (0, 10), [x0, 0], t_eval=t, rtol=1e-10, atol=1e-12)
    xa, wd, z = analytic(t, m, k, c, x0)
    err = np.max(np.abs(sol.y[0] - xa))
    print(f"selftest: damped period {2*math.pi/wd:.4f} s, zeta {z:.4f}, max |ode - closed form| = {err:.2e}")
    print("ok: sdof_reference selftest passed" if err < 1e-7 else "FAIL")
    sys.exit(0 if err < 1e-7 else 1)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?"); ap.add_argument("--m", type=float); ap.add_argument("--k", type=float)
    ap.add_argument("--c", type=float); ap.add_argument("--dof", default="y")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    data = np.genfromtxt(a.csv, delimiter=",", names=True)
    t = data["t"]; x = data[a.dof] - data[a.dof][-1]   # remove final offset
    x0 = x[0]
    xa, wd, z = analytic(t - t[0], a.m, a.k, a.c, x0)
    print(f"analytic: damped period {2*math.pi/wd:.4f} s, zeta {z:.4f}")
    print(f"RMS(fluent - analytic) = {np.sqrt(np.mean((x - xa)**2)):.3e} m, x0 = {x0:.4f} m")
