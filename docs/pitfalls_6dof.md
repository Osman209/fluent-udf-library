# Six-DOF pitfalls (from runs, not from the manual)

These are the things that cost days on a floating-hull VOF six-DOF case
in Fluent 2025 R1. They are not in the UDF manual in this form.

## 1. Only one zone may be the active six-DOF body

If you use an overset mesh, the hull wall **and** the overset component
zone both need Six DOF motion, but only the hull wall is the body. The
overset zone must be set to **Six DOF + Passive**. With both set as
active bodies Fluent prints

    Info: 6DOF: can't compute angular acceleration. Check moments/products of inertia.

on every iteration, the rotational DOFs are dead, and the body
translates but never turns. This message is a setup error, not an
inertia value error.

## 2. The mass must be the real displacement mass

A placeholder mass (say 1000 kg for a hull that displaces 8000 kg) makes
the run diverge in the first second with the mesh-motion residual
growing inside the time step (28 -> 75 -> 250 in three iterations) and
continuity stuck. This is the added-mass instability: the fluid force
on the body is proportional to its acceleration, and if the body is
much lighter than the water it displaces, the explicit coupling is
unstable. Use the real mass and inertia. If you have to run a light
body, reduce the time step and increase the six-DOF sub-iterations.

## 3. Let the hull settle before you push it

Start with a float phase: no thrust, no turn, body free in heave, roll
and pitch, for at least two natural heave periods (for an 11 m hull that
is 3 to 5 seconds). The controller UDF has `T_HOLD` for this. If heave
drifts far from the reference draft in the float phase (we saw +36 cm),
the mass or the centre of gravity is not the real one; fix that before
going on.

## 4. Check the frozen-body case first

Run the same case with the body frozen (all DOFs locked). If it does
not converge cleanly, the problem is the mesh, the overset setup or the
boundary conditions, not the six-DOF coupling. If the frozen case is
clean and the free case diverges, it is one of items 1 to 3.

## 5. Macro name changed: DT_OMEGA

In Fluent 2025 R1, `DT_OMEGA(dt)` is a **scalar** (angular speed about
the rotation axis of a rigid rotation). `DT_OMEGA(dt)[1]` does not
compile. The angular velocity **vector** of a six-DOF body is
`DT_OMEGA_CG(dt)`. `DT_CG`, `DT_VEL_CG` and `DT_THETA` are arrays as
expected. This library uses `DT_OMEGA_CG` everywhere.

## 6. The SDOF UDF is called more than once per time step

`DEFINE_SDOF_PROPERTIES` runs inside the six-DOF sub-iterations. If you
log from it, guard with "only when the flow time has advanced" (see
`sdof_spring_damper.c`), otherwise you get several lines per step and a
file that is hard to plot.

## 7. Gravity is already in

If gravity is on in Operating Conditions, the six-DOF solver applies
the weight. Do not add m g in the UDF as well; the body will sink at
twice the rate and you will blame the mesh.

## 8. Reversed flow at the outlet is a symptom, not the cause

"Reversed flow in N faces on pressure-outlet" appearing together with a
diverging six-DOF run is almost always caused by the body motion, not
by the outlet. Fix the motion first.
