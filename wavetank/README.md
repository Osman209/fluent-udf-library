# wavetank pack

| File | Macro | Hook |
|---|---|---|
| wave_inlet_regular.c | DEFINE_PROFILE x3 | velocity inlet: X-Vel, Y-Vel, water Volume Fraction |
| wave_inlet_irregular.c | DEFINE_PROFILE x3 | same |
| damping_zone_source.c | DEFINE_SOURCE | Cell Zone > Source Terms > momentum (mixture). In 2D hook x and y only |
| wave_probe.c | DEFINE_EXECUTE_AT_END | Function Hooks > Execute at End |
| wave_init.c | DEFINE_INIT | Function Hooks > Initialization |

## Tank setup that goes with these files

- 2D or 3D, x along the tank, y up, still water level at Y_SWL (the same
  value in every file), bed at Y_SWL - DEPTH_D.
- Multiphase: VOF, implicit body force on, air primary phase, water
  secondary (so WATER_PHASE = 1). Operating density = air density.
- Inlet: velocity inlet on the whole left wall (air and water), with the
  three profiles hooked. Top: pressure outlet.
- Far end: a **wall**, with the damping zone in front of it. A pressure
  outlet there does not balance the hydrostatic column of water, so the
  tank drains: with a 1 m depth and the operating density set to air,
  the gauge pressure at the bed is about 9.8 kPa pushing against 0 Pa.
  The beach absorbs the wave before it reaches the wall, so the wall
  reflects almost nothing.
- Mesh: at least 10 cells per wave height vertically near the free
  surface, 60 to 100 cells per wavelength horizontally. Time step small
  enough that the Courant number at the surface stays below about 0.25.
- Length: about 3 wavelengths of working section plus 1.5 to 2
  wavelengths of beach (X_START to X_END in damping_zone_source.c).

## Validating the tank (once per mesh and time step)

0. Before anything else: run one time step and look at `wave_probes.csv`.
   On still water every probe must read close to zero. A constant offset
   means the probe geometry is wrong, and every later number inherits it.
1. Run with no body, regular wave, 15 to 20 periods.
2. `python3 validation/wave_probe_analysis.py wave_probes.csv --H 0.10 --T 1.5 --d 1.0`
3. The wave height at the first probe should be within about 5 % of the
   target after the ramp, and the period within about 1 %. The reflection
   coefficient from two probes in the working section should be low; if
   it is not, lengthen the beach or raise C0.
   The first two probes are the pair used for reflection, and they must
   be between 0.05 and 0.45 wavelengths apart. Change the wave and you
   must move them: the script refuses to report a coefficient outside
   that range rather than print a wrong one.
4. Then put the body in.

The regular inlet prints the Ursell number at start-up. Above about 26,
Stokes second order is outside its range.
