"""
check_defaults.py

The values shipped in the files are what someone gets if they build the
library without editing anything. They have to describe one coherent,
physically valid case, or the first run teaches the wrong lesson.

This caught a probe pair shipped 0.90 wavelengths apart: outside the
0.05-0.45 window the two-probe reflection method needs, and far enough
from a multiple of pi that the script's own sin(k dx) warning stayed
quiet. The reflection coefficient came out wrong and looked fine.
"""
import re, math, sys, os

ROOT = os.path.join(os.path.dirname(__file__), "..")

def defs(rel):
    d = {}
    for m in re.finditer(r'#define\s+([A-Z][A-Z0-9_]*)\s+(-?[\d.]+(?:e-?\d+)?)\s',
                         open(os.path.join(ROOT, rel), encoding='utf-8').read()):
        try: d[m.group(1)] = float(m.group(2))
        except ValueError: pass
    return d

def wavelength(T, d, g=9.81):
    w = 2 * math.pi / T
    k = w * w / g
    for _ in range(200):
        k -= (g * k * math.tanh(k * d) - w * w) / (g * math.tanh(k * d) + g * k * d / math.cosh(k * d) ** 2)
    return 2 * math.pi / k

bad = []
def need(c, msg):
    if not c: bad.append(msg)

reg = defs('wavetank/wave_inlet_regular.c')
dmp = defs('wavetank/damping_zone_source.c')
prb = defs('wavetank/wave_probe.c')
ini = defs('wavetank/wave_init.c')
irr = defs('wavetank/wave_inlet_irregular.c')

need(reg['Y_SWL'] == prb['Y_SWL'] == ini['Y_SWL'], "still water level differs between files")
need(prb['Y_BED'] == reg['Y_SWL'] - reg['DEPTH_D'], "probe bed level does not match SWL minus depth")
need(prb['Y_TOP'] > prb['Y_SWL'], "probe domain top is not above the water line")
need(reg['DEPTH_D'] == irr['DEPTH_D'], "regular and irregular inlets assume different depths")

L = wavelength(reg['WAVE_T'], reg['DEPTH_D'], reg['GRAV'])
H = reg['WAVE_H']
need(H * L * L / reg['DEPTH_D'] ** 3 < 26, "the default wave is outside the Stokes second-order range")
need(H / L < 0.14, "the default wave is steeper than the breaking limit")

beach = (dmp['X_END'] - dmp['X_START']) / L
need(1.4 < beach < 2.2, f"the default beach is {beach:.2f} wavelengths, not the 1.5-2 the file advises")

px = [float(x) for x in re.search(r'PROBE_X\s+\{([^}]+)\}',
      open(os.path.join(ROOT, 'wavetank/wave_probe.c'), encoding='utf-8').read()).group(1).split(',')]
need(all(x < dmp['X_START'] for x in px), "a default probe sits inside the beach")
gap = (px[1] - px[0]) / L
need(0.05 < gap < 0.45,
     f"the first probe pair is {gap:.2f} wavelengths apart; the two-probe reflection "
     f"method needs 0.05 to 0.45, and outside it the answer is wrong without looking wrong")

need(reg['VOF_SMOOTH'] < H / 2, "the inlet smoothing band is not small against the wave height")
need(ini['SMOOTH_DY'] < H / 2, "the initialisation smoothing is not small against the wave height")
need(reg['RAMP_PERIODS'] >= 1, "the inlet ramp is shorter than one period")
need(irr['RAMP_S'] >= reg['WAVE_T'], "the irregular inlet ramp is shorter than one period")

if bad:
    print("FAIL: the shipped defaults do not describe a coherent case")
    for b in bad: print("  " + b)
    sys.exit(1)
print(f"ok: the shipped defaults describe one coherent tank "
      f"(L = {L:.2f} m, beach {beach:.2f} L, probe pair {gap:.2f} L)")
