# Unreleased - core validation fixes

- Correct piecewise-linear cp enthalpy integration and boundary continuation.
- Correct fluid-on-wall viscous force sign in gradient method.
- Match restart state pointers in single/double precision; exclude host writer.
- Strict numeric 1D tables with finite/increasing checks and row diagnostics.
- Allocation and 2D table finite/grid-size checks.
- Explicit fixed spring datum; positive-time spring logs append on restart.
- Actual-UDF regressions, 204 mock syntax configurations, CI, catalog and
  a wave-tank acceptance protocol. No new Fluent runtime validation claimed.

# Changelog

## 0.1.0 - 2026-09-10

First release. Three packs: motion (7 UDFs), wave tank (5 UDFs) and
profiles (5 UDFs). Shared headers for parallel-safe file writing, table
reading and wave theory. Test suite that syntax-checks every UDF in 2D
and 3D and cross-checks the physics against independent numpy
derivations. Docs on compiling, parallel safety and six-DOF pitfalls.
