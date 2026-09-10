# motion pack

| File | Macro | Hook |
|---|---|---|
| cg_motion_translate.c | DEFINE_CG_MOTION | Rigid Body > Motion UDF |
| cg_motion_oscillate.c | DEFINE_CG_MOTION | Rigid Body > Motion UDF |
| cg_motion_from_table.c | DEFINE_CG_MOTION, DEFINE_ON_DEMAND | Rigid Body > Motion UDF |
| sdof_spring_damper.c | DEFINE_SDOF_PROPERTIES | Rigid Body > Six DOF UDF |
| sdof_speed_heading_controller.c | DEFINE_SDOF_PROPERTIES | Rigid Body > Six DOF UDF |
| grid_motion_flapping_plate.c | DEFINE_GRID_MOTION | User-Defined > Mesh Motion UDF |
| force_moment_logger.c | DEFINE_EXECUTE_AT_END, DEFINE_ON_DEMAND | Function Hooks > Execute at End |

## Which one do I need

- Body moves how I tell it (forced motion): `cg_motion_*`. Fluent moves
  the mesh; the fluid does not push back on the body.
- Body moves because the fluid pushes it (free motion, VIV, floating
  hull): `sdof_*`. Fluent integrates the equations of motion with the
  fluid loads plus what you add in the UDF.
- Surface bends (plate, membrane): `grid_motion_*`.
- Forces in a file every time step, in parallel: `force_moment_logger.c`.

## Checking a run

- Forced motion: read the body position back from Fluent (Report >
  Dynamic Mesh, or the CG in the console) and compare with the law in the
  file. `validation/motion_reference.py` gives the exact position laws.
- Free motion: run a free-decay test (release the body from an offset,
  fluid loads small) and compare `sdof_motion.csv` with
  `validation/sdof_reference.py`. The damped period should match
  2 pi / (wn sqrt(1 - zeta^2)) before you trust the coupled run.
- Forces: run `calibrate_viscous_sign` from Execute on Demand once on a
  converged solution and compare with Report > Forces on the same zone.
  It prints the pressure force and the viscous force computed two ways.

Before anything else on a six-DOF VOF case, read
`../docs/pitfalls_6dof.md`.
