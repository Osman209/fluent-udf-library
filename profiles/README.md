# profiles pack

Boundary conditions and material properties.

| File | Macro | Hook |
|---|---|---|
| profile_abl_wind.c | DEFINE_PROFILE x5 | velocity inlet: velocity components, k, epsilon or omega |
| profile_parabolic.c | DEFINE_PROFILE x3 | velocity inlet: normal velocity component |
| profile_transient_csv.c | DEFINE_PROFILE x3, DEFINE_ON_DEMAND | any boundary quantity set to user-defined |
| profile_pulsatile.c | DEFINE_PROFILE x2 | velocity inlet |
| property_from_table.c | DEFINE_PROPERTY x4, DEFINE_SPECIFIC_HEAT | Materials panel, property set to user-defined |

## Notes worth knowing before you use them

**ABL wind.** The k and epsilon profiles are the Richards-Hoxey set that
is consistent with the log law, so the inlet profile is a solution of the
turbulence model and does not decay down the domain by itself. Keeping it
from decaying also needs a ground roughness that matches z0 and a top
boundary that does not fight the profile; that is case setup, not UDF.
Roughness by terrain: 0.0002 m open sea, 0.03 m open farmland, 0.3 m
suburban, 1.0 m dense urban.

**Parabolic and power-law inlets.** All three are normalised so the
area-average velocity is exactly the U_MEAN you set, whatever the mesh.
So the check is simple: the inlet mass flow should equal rho * U_MEAN *
inlet area. If it does not, the profile is hooked to the wrong component
or the cross-section indices are wrong.

**Transient from file.** DEFINE_PROFILE is re-evaluated every iteration,
so the value tracks the flow time in a transient run. In a steady run
the flow time does not advance and the value stays at the first row.
Interpolation is linear between the rows you give, so give enough rows
for the shape you want.

**Pulsatile.** A parabolic profile is only right at low Womersley number
(below about 1). For an aorta, alpha is roughly 10 to 20 and the real
profile is much flatter. The UDF prints alpha at start-up and warns when
the flat version is the better choice.

**Specific heat.** DEFINE_SPECIFIC_HEAT must return cp AND the sensible
enthalpy measured from the reference temperature; Fluent uses that
enthalpy in the energy equation. This file integrates your cp table once
at load time so the two cannot disagree. Supplying an enthalpy that does
not match the cp is the usual reason a cp UDF gives a wrong temperature
field while cp itself looks right.

## File formats

1D property or transient table, two columns, `#` for comments, comma,
space or tab separated:

    # T [K]   viscosity [Pa s]
    273.15    0.001792
    283.15    0.001307

2D property table (temperature and pressure), a regular grid:

    nT nP
    T1 T2 ... TnT
    P1 P2 ... PnP
    v(T1,P1) v(T1,P2) ... v(T1,PnP)
    v(T2,P1) ...
