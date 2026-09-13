# Validation status and release gate

This revision was checked outside Fluent. No new Fluent build, serial run,
parallel run or experimental validation is claimed. Original notes mention
Fluent 2025 R1; they are historical observations, not certification of this
revision. See `catalog/udfs.json` for the inventory of all 17 source files.

## Evidence actually obtained

- 204 mock syntax configurations: 17 files x 2 dimensions x 2 precisions x
  3 process roles (serial, host, compute node).
- Actual `cp_T` execution: off-knot enthalpy, dh/dT, both extrapolation
  directions, variable/reversed reference, and single-point table.
- Actual `visc_force_gradient`: a controlled bottom-wall shear field and
  its reverse; checked in 2D/3D and float/double mocks.
- Actual `spring_damper`: displaced first call, time reset and simulated
  library reload retain the fixed physical datum.
- Strict reader rejects text headers, missing/extra fields, empty CSV fields,
  nonfinite values, overflow, and non-increasing axes.
- Writer selection: serial=yes, host=no, node 0=yes, other node=no.
- Existing dispersion, wave, spectrum, motion and profile reference checks.

Mock tests do not reproduce MPI, solver memory allocation, macro ABI,
mesh updates, VOF evolution, or wall-function physics. They do not prove
that every listed configuration is supported in Fluent.

## Required Fluent evidence before a validated release

For each UDF record the exact version/build, OS/compiler, dimension, precision,
solver and models, hook settings, serial/parallel process count, input files,
mesh, time step, reference and acceptance tolerance. Attach console and results.
Use `examples/wavetank/README.md` for the first end-to-end validation pack.

Proposed labels: `mock-checked`, `fluent-build-verified`,
`fluent-case-verified`. Advance only with archived evidence. Formula-only
checks must not be presented as tests of a different implementation.

## Known limitations remaining

- Force logger method 1 needs available velocity gradients and comparison
  with Fluent wall forces. Method 2 uses undocumented wall-shear storage;
  its calibration hook still references that storage even with method 1.
- Wave probes assume a rectangular constant-height water column; variable
  bathymetry, excluded fluid volume or solids invalidate that conversion.
  Empty probes currently return zero; inspect column coverage separately.
- Table readers assume all compute processes can access the same input files.
  Missing/rejected tables retain the documented legacy fallback properties
  or zero motion/BC, with a diagnostic; inspect it before running a case.
- The speed/heading controller uses a simplified yaw-only forward vector;
  large roll/pitch and wrapped headings need a separate controller update.
- Axis choices and mesh quadrature still need checks. An analytically
  normalized inlet is not guaranteed to have exactly normalized discrete flux.
- Clock reversal is a logging heuristic, not a complete restart detector.
  Most loggers still overwrite on library reload. Spring logs append at a
  positive restart time and can repeat the saved-time row; separate run
  directories are recommended.
