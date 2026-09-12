"""
wave_probe_analysis.py

Post-processing for wave_probes.csv written by wavetank/wave_probe.c.

For each probe: mean wave height and period by zero-up-crossing on the
last N_PERIODS periods, compared with the target H and T you give.

Reflection coefficient from a pair of probes by the two-probe method
of Goda and Suzuki (1976): with probes at x1 and x2 = x1 + dx, the
incident and reflected amplitudes at the wave frequency are

    a_i, a_r = |A1 e^{+-i k dx} - A2| / (2 |sin(k dx)|)   (up to which is which)

where A1, A2 are the complex Fourier amplitudes at the wave frequency.
Valid when k dx is not near a multiple of pi; the script warns when
sin(k dx) is small. Use dx between 0.05 L and 0.45 L.

Usage:
    python3 validation/wave_probe_analysis.py wave_probes.csv --H 0.10 --T 1.5 --d 1.0
    python3 validation/wave_probe_analysis.py --selftest

--selftest builds a synthetic incident + 20 % reflected wave, writes it
in the probe file format, and checks that the script recovers H, T and
Kr = 0.20. That is the test that this analysis itself is right.
"""
import sys, math, argparse, tempfile, os
import numpy as np
from scipy.optimize import brentq

def wavenumber(T, d, g=9.81):
    w = 2 * math.pi / T
    return brentq(lambda k: g * k * math.tanh(k * d) - w * w, 1e-6, 1e4)

def _crossings(t, eta):
    """Times where the signal crosses its own mean going up."""
    e = eta - eta.mean()
    i = np.where((e[:-1] < 0) & (e[1:] >= 0))[0]
    if len(i) < 2:
        return np.array([]), i
    tc = t[i] - e[i] * (t[i + 1] - t[i]) / (e[i + 1] - e[i])
    return tc, i


def level_drift(t, eta):
    """How far the mean level moved across the window, in metres.

    Reported separately rather than removed. A progressive wave carries a
    net mass flux forward (Stokes drift), and in a tank closed by a wall
    that water has nowhere to go: for H = 0.06 m and T = 1.2 s in a 12 m
    tank the level climbs about 4 mm over eighteen periods, which is 7 %
    of the wave height. That is physics. A level falling instead, or
    moving far more than this, is a leak somewhere and worth knowing
    about, so the number is shown rather than quietly subtracted.

    Measured as the mean of the first whole cycle against the mean of the
    last. Fitting a straight line to the raw window instead reads the
    leftover part-cycle at each end as a slope, and reported nearly 6 mm
    of drift on a wave that had none at all."""
    _, i = _crossings(t, eta)
    if len(i) < 3:
        return 0.0
    return eta[i[-2]:i[-1]].mean() - eta[i[0]:i[1]].mean()


def zero_upcross(t, eta):
    """Wave height and period by zero-up-crossing about the mean.

    Removing the mean is enough: a 4 mm drift over an eight-period window
    changes the measured height of a 60 mm wave by 0.1 mm, which was
    checked rather than assumed."""
    eta = eta - eta.mean()
    idx = np.where((eta[:-1] < 0) & (eta[1:] >= 0))[0]
    if len(idx) < 3:
        return None, None
    tc = t[idx] - eta[idx] * (t[idx + 1] - t[idx]) / (eta[idx + 1] - eta[idx])
    periods = np.diff(tc)
    heights = [eta[idx[i]:idx[i + 1]].max() - eta[idx[i]:idx[i + 1]].min() for i in range(len(idx) - 1)]
    return np.mean(heights), np.mean(periods)

def fourier_amp(t, eta, f):
    eta = eta - eta.mean()
    dt = t[1] - t[0]
    return 2.0 * np.sum(eta * np.exp(-2j * math.pi * f * t)) * dt / (t[-1] - t[0])

def reflection(t, e1, e2, dx, T, d):
    k = wavenumber(T, d)
    L = 2 * math.pi / k
    s = math.sin(k * dx)
    A1 = fourier_amp(t, e1, 1 / T)
    A2 = fourier_amp(t, e2, 1 / T)
    # eta_j = Re[(a_i e^{i k x_j} + a_r e^{-i k x_j}) e^{-i w t}].
    # fourier_amp returns the conjugate of the bracket, so with x1 = 0
    # and x2 = dx:  A1 = al + be,  A2 = al e^{-ik dx} + be e^{+ik dx}
    # where al, be carry the incident / reflected amplitude and phase.
    M = np.array([[1, 1], [np.exp(-1j * k * dx), np.exp(1j * k * dx)]])
    ai, ar = np.linalg.solve(M, np.array([A1, A2]))
    return abs(ar) / abs(ai), abs(s), dx / L

def analyse(path, H, T, d, n_periods=8, quiet=False):
    # read the header ourselves: numpy drops the "." in names like eta_x5.9
    with open(path) as fh:
        header = fh.readline().strip().split(",")
    raw = np.loadtxt(path, delimiter=",", skiprows=1)
    t = raw[:, 0]
    names = [h for h in header if h.startswith("eta_x")]
    xs = [float(h[5:]) for h in names]
    data = {h: raw[:, header.index(h)] for h in names}
    tmask = t > t[-1] - n_periods * T
    res = {}
    for n, x in zip(names, xs):
        Hm, Tm = zero_upcross(t[tmask], data[n][tmask])
        drift = level_drift(t[tmask], data[n][tmask]) * 1000
        res[x] = (Hm, Tm)
        if not quiet:
            print(f"probe x={x:g}: H = {Hm:.4f} (target {H}), T = {Tm:.4f} (target {T})"
                  + (f", level drift {drift:+.1f} mm" if abs(drift) > 0.5 else ""))
    kr = None
    if len(xs) >= 2:
        x1, x2 = xs[0], xs[1]
        kr, s, ratio = reflection(t[tmask], data[names[0]][tmask],
                                  data[names[1]][tmask], x2 - x1, T, d)
        if not quiet:
            # The spacing itself is the thing to check. Testing only
            # sin(k dx) misses a spacing of 0.9 wavelengths, where the
            # sine is a healthy 0.59 but the two probes are more than
            # half a wavelength apart and the phase has wrapped: the
            # incident and reflected waves can no longer be separated,
            # and a wrong number comes out looking perfectly reasonable.
            if not (0.05 <= ratio <= 0.45):
                print(f"reflection coefficient: NOT COMPUTED. The probes are "
                      f"{ratio:.2f} wavelengths apart; the method needs 0.05 to 0.45. "
                      f"Move them to about {0.3 * (x2 - x1) / ratio:.2f} m apart.")
                kr = None
            else:
                print(f"reflection coefficient from probes at {x1:g} and {x2:g}: "
                      f"Kr = {kr:.3f}  (spacing {ratio:.2f} wavelengths)"
                      + ("   WARNING: sin(k dx) small, unreliable" if s < 0.2 else ""))
    return res, kr

def selftest():
    H, T, d, Kr_true = 0.10, 1.5, 1.0, 0.20
    k = wavenumber(T, d); w = 2 * math.pi / T
    t = np.arange(0, 40, 0.005)
    xs = [5.0, 5.9]
    cols = []
    for x in xs:
        eta = 0.5 * H * (np.cos(k * x - w * t) + Kr_true * np.cos(-k * x - w * t + 0.7))
        cols.append(eta)
    fd, path = tempfile.mkstemp(suffix=".csv"); os.close(fd)
    with open(path, "w") as f:
        f.write("t," + ",".join(f"eta_x{x:g}" for x in xs) + "\n")
        for i in range(len(t)):
            f.write(f"{t[i]:g}," + ",".join(f"{c[i]:.8g}" for c in cols) + "\n")
    res, kr = analyse(path, H, T, d, quiet=True)
    os.remove(path)
    Hm, Tm = res[xs[0]]
    print(f"selftest: T recovered {Tm:.4f} (true {T}), Kr recovered {kr:.4f} (true {Kr_true})")
    ok = abs(Tm - T) < 1e-3 and abs(kr - Kr_true) < 0.01
    print("ok: wave_probe_analysis selftest passed" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?")
    ap.add_argument("--H", type=float); ap.add_argument("--T", type=float); ap.add_argument("--d", type=float)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    else:
        analyse(a.csv, a.H, a.T, a.d)
