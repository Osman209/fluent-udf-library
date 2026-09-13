# fluent-udf-library

User Defined Function (UDF) examples for ANSYS Fluent, in three packs.

**Status:** standalone checks pass; this revision has not been built or run
in Fluent here. See [validation status](docs/validation-status.md) and the
[machine-readable catalog](catalog/udfs.json) before choosing a UDF.

**motion/** - rigid-body and deforming-mesh motion: prescribed translation
and oscillation, motion from a table, six-DOF with springs and dampers, a
speed-and-heading controller for a self-propelled hull, a deforming
flapping plate, and a force/moment logger.

**wavetank/** - a VOF numerical wave tank: regular waves (Airy and Stokes
second order), irregular waves from a JONSWAP spectrum, a damping zone
(numerical beach), wave probes, and still-water initialisation.

**profiles/** - boundary conditions and material properties: atmospheric
boundary layer wind with matching turbulence, fully developed pipe and
channel inlets, any boundary condition varying with time from a file,
pulsatile inlets, and material properties interpolated from data.

Every file has its user parameters at the top, a comment saying what it
does and where to hook it. Serial and parallel use must be verified for
your Fluent version and case; mock compilation is not solver validation.

## Why this exists

Most UDF code online is a single file someone pasted once. The usual
problems: it only works in serial, it double counts faces in parallel, it
writes the log file from every process, it starts the motion with an
impulse and the run blows up on the first step, or a formula has a sign
wrong. This library handles those once and keeps a test for each.

[ملاحظات الترقية بالعربية](docs/UPGRADE-AR.md)

## Quick start

1. Copy the `.c` file you need, plus the headers from `common/` that it
   includes, into your Fluent working directory.
2. Edit the parameters at the top of the `.c` file.
3. Fluent: User-Defined > Functions > Compiled > add the source file (and
   the headers under "Header Files") > Build > Load.
4. Hook it where the file's header comment says.

Use **compiled**, not interpreted: several files use `dynamesh_tools.h`,
file I/O and static arrays, which the interpreter does not support.
`docs/compile_checklist.md` covers the build errors people hit most.

## Tests

    sh tests/run_tests.sh

Needs gcc and Python 3. Install reference dependencies with
`python3 -m pip install -r validation/requirements.txt`. It runs:

- a gcc mock syntax matrix for 2D/3D, float/double, serial/host/node
- regressions that call the shipped cp, shear-force and spring UDF code
- strict table validation and writer-selection tests
- plain-C tests of the shared wave theory, the table reader and the
  profile formulas
- numpy cross-checks that re-derive the same physics independently:
  Airy and Stokes second order agree to 1e-16, the JONSWAP components
  match to 1e-13, the ABL profiles satisfy the Richards-Hoxey
  consistency condition, the pipe and channel profiles integrate to the
  mean velocity you asked for, the cp table gives an enthalpy consistent
  with itself, and the two-probe reflection analysis recovers a known
  reflection coefficient from a synthetic signal

## Files

    common/
      udf_common.h          parallel-safe writer macro, interpolation, table reader
      wave_theory.h         dispersion, Airy, Stokes 2nd, JONSWAP, irregular sea
    motion/
      cg_motion_translate.c            constant velocity with smooth ramp
      cg_motion_oscillate.c            harmonic heave and pitch
      cg_motion_from_table.c           motion from a time table file
      sdof_spring_damper.c             6DOF with springs, dampers, locks, motion log
      sdof_speed_heading_controller.c  thrust + yaw controller for a hull
      grid_motion_flapping_plate.c     deforming plate, parallel-safe node update
      force_moment_logger.c            forces and moments on a wall
    wavetank/
      wave_inlet_regular.c    Airy / Stokes 2nd order inlet (u, v, VOF)
      wave_inlet_irregular.c  JONSWAP irregular inlet
      damping_zone_source.c   momentum sink beach
      wave_probe.c            free-surface elevation at stations
      wave_init.c             still-water VOF initialisation
    profiles/
      profile_abl_wind.c      ABL log / power law + k, epsilon, omega
      profile_parabolic.c     laminar pipe, laminar channel, 1/7 power law pipe
      profile_transient_csv.c any BC varying with time, from a file
      profile_pulsatile.c     Fourier-series pulsatile inlet
      property_from_table.c   rho, mu, k, cp from data; cp with consistent enthalpy
    validation/               reference solutions and post-processing scripts
    tests/                    gcc tests and the mock udf.h
    docs/                     compiling, parallel safety, six-DOF pitfalls

## If something does not work

Fluent changes macros between versions, and every case is set up a little
differently. If a file does not build, does not hook, or gives a result
you do not expect, get in touch and it will be sorted out:

- Open an issue on GitHub. There are templates for a build problem, an
  unexpected result, and a request for a UDF that is not here yet; they
  ask for the Fluent version, the operating system and the console
  output, which is usually enough to find the cause first time.
- Or email mohamedosmannn2999@gmail.com if you would rather not post
  publicly.

Pull requests are welcome, including new UDFs for other physics. See
`CONTRIBUTING.md` for the house style and the checks.

## Licence

MIT. Use it in your work, commercial or not. Attribution appreciated, not
required. As with any CFD setup, check the results against a case you know
before you rely on them.

## Author

Osman - mechanical engineer, CFD with ANSYS Fluent, STAR-CCM+ and LS-DYNA.
Arabic CFD tutorials on YouTube. mohamedosmannn2999@gmail.com

## Upgrade notes: corrected core

- Set `SPRING_REF_X` and `SPRING_REF_TH` to the physical equilibrium.
  The spring no longer captures the first observed position. Defaults are zero.
- Numeric table headers must start with `#`. Bad rows now reject the entire
  table, rather than being skipped or converted to zero. First column must
  strictly increase; values must be finite. Legacy fallback behavior remains
  and is printed; verify table loading before starting a case.
- Specific heat now integrates the interpolated cp exactly between knots
  and continues the enthalpy outside the table with clamped endpoint cp.
- Default viscous force now reports fluid-on-wall direction. Recheck
  force and moment comparisons before reusing old calibration results.
- Begin Fluent verification with [the reference tank](examples/wavetank/README.md).
