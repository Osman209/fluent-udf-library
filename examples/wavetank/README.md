# First reference tank: acceptance protocol (not yet Fluent-verified)

A reproducible setup specification, not a prebuilt .cas/.dat or tested
version-specific journal. Record the actual Fluent settings and mesh used.

## Baseline configuration

- 2D planar rectangular tank: x=0..25 m, y=-1..0.5 m; y up.
- Air primary, water secondary (index 1), VOF, gravity (0,-9.81) m/s2.
- Regular wave: H=0.10 m, T=1.50 s, d=1 m, still-water y=0;
  Stokes second order, ramp=2 periods; g=9.81 m/s2.
- Use velocity components at x=0: `inlet_u`, `inlet_v`; water VOF:
  `inlet_vof`. Initialize with `still_water` after setting material models.
- Bottom: wall. Top and outlet: atmospheric pressure with appropriate
  phase backflow; verify hydrostatic outlet treatment for the selected
  Fluent version. Do not assume this specification alone prevents outlet artifacts.
- Damping source: x=20..25 m, C0=10 1/s, exponent=2; mixture x/y momentum.
- Probes: x=5,6,10,14 m; width=0.10 m, y bed=-1, top=0.5, water index=1.
- Reference dispersion: k=1.87477235386 1/m, L=3.35143906631 m.

Start with a resolved free-surface mesh and at least 100 time steps per
period; these are starting choices, not convergence evidence. Refine both.

## Run sequence and proposed gates

1. Still water with wave amplitude zero: verify hydrostatic balance and
   probe error below 0.25 local vertical cell heights after settling.
2. Regular waves: sample after ramp AND propagation to the gauges. Retain
   at least 10 clean periods and exclude initial transients/reflection arrival
   as appropriate for the measurement. Compare period (target <=1% error),
   wave height (<=5%) and phase against the supplied theory.
3. Use the first probe pair for reflection analysis. The spacing is 0.30 L;
   proposed beach target Kr<=0.10. Tune and repeat if it fails; do not change
   the target retroactively to label a case as passing.
4. Halve time step, then refine the interface/propagation mesh: seek <=2%
   change in the chosen gauge amplitude and report phase/reflection changes.
5. Repeat the same case in serial, 2 and 4 compute processes. Compare
   phase-aligned signals and integral mass balance with stated tolerances;
   tiny floating-point differences alone do not prove a parallel bug.
6. Save at a displaced surface and resume. Compare with an uninterrupted run;
   keep log outputs separate to avoid overwriting or duplicate sample times.

These gates are proposed acceptance criteria, not measured results.
Use `validation/wave_probe_analysis.py --help` for the existing analysis CLI.

## Evidence record

Copy `validation-record.json` per run and fill actual values. Archive mesh,
settings, hook names, console transcript, input revision and CSV outputs.
Never add a pass badge without those files or an accessible evidence link.
