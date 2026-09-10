"""
profiles_check.py

Independent numpy check of the profile formulas in profiles/.
Re-derives each one from its definition rather than copying the code:

  - ABL log law: fits U(z) over the domain and recovers ustar from the
    slope against ln(z+z0); checks it equals kappa U_REF / ln((Zref+z0)/z0)
  - ABL turbulence: checks d(nu_t)/dz = kappa ustar, the condition for
    the log law to be a solution of the k-epsilon equations
  - pipe and channel profiles: area-averages by quadrature
  - power-law pipe: checks the closed-form ratio
    U_mean/U_max = 2n^2/((n+1)(2n+1)) against direct integration

Run:  python3 validation/profiles_check.py
"""
import numpy as np
from scipy.integrate import quad

kappa, Cmu, Uref, Zref, z0 = 0.42, 0.09, 10.0, 10.0, 0.03
ustar = kappa * Uref / np.log((Zref + z0) / z0)

z = np.linspace(0.01, 200, 5000)
U = (ustar / kappa) * np.log((z + z0) / z0)
slope = np.polyfit(np.log(z + z0), U, 1)[0]
print(f"ABL: ustar from the slope = {slope*kappa:.8f}, direct = {ustar:.8f}")
assert abs(slope * kappa - ustar) / ustar < 1e-10

k = ustar**2 / np.sqrt(Cmu)
eps = ustar**3 / (kappa * (z + z0))
nut = Cmu * k**2 / eps
dnut = np.gradient(nut, z)
print(f"ABL: d(nu_t)/dz = {dnut.mean():.8f}, kappa*ustar = {kappa*ustar:.8f}")
assert abs(dnut.mean() - kappa * ustar) / (kappa * ustar) < 1e-6

R, Um, n = 0.05, 1.0, 7.0
num = quad(lambda r: 2 * Um * (1 - r**2 / R**2) * 2 * np.pi * r, 0, R)[0]
print(f"parabolic pipe area-average = {num / (np.pi * R**2):.10f}")
assert abs(num / (np.pi * R**2) - Um) < 1e-10

umax = Um * ((n + 1) * (2 * n + 1)) / (2 * n**2)
num = quad(lambda r: umax * (1 - r / R) ** (1 / n) * 2 * np.pi * r, 0, R)[0]
print(f"power-law pipe area-average = {num / (np.pi * R**2):.10f} (umax = {umax:.6f})")
assert abs(num / (np.pi * R**2) - Um) < 1e-8

H = 0.05
num = quad(lambda y: 1.5 * Um * (1 - y**2 / H**2), -H, H)[0]
print(f"parabolic channel average = {num / (2 * H):.10f}")
assert abs(num / (2 * H) - Um) < 1e-10

print("ok: all profile formulas match independent numpy derivations")
