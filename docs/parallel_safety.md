# Parallel safety

A parallel Fluent run has one **host** process (no mesh, does the GUI and
file I/O) and N **compute nodes** (each has a partition of the mesh).
Serial Fluent is a single process that behaves like one compute node.
The same UDF source is compiled for all of them, and the macros `RP_HOST`
and `RP_NODE` tell you which one is running the code.

Three rules cover almost every parallel failure:

## 1. Do not loop over the mesh on the host

    #if !RP_HOST
       ... begin_c_loop / begin_f_loop ...
    #endif

The host has no cells; a loop there reads invalid memory and you get a
segmentation fault that does not happen in serial. Every loop in this
library is wrapped this way. `DEFINE_PROFILE` is an exception: Fluent
only calls it on the nodes.

## 2. Sum, then reduce, then let one process write

    Fx = PRF_GRSUM1(Fx);          /* global sum across nodes */
    if (UDF_IS_WRITER) fprintf(...)

`PRF_GRSUM1` adds the partial sums from all nodes. Without it each node
logs only its own partition. `UDF_IS_WRITER` (in `udf_common.h`) is false on the host, true
on compute node 0 in parallel and true in serial, so the file gets one
line per step instead of N.

For faces on partition boundaries, a face can exist on two nodes. Use
`PRINCIPAL_FACE_P(f,t)` to count it once. For cells use
`begin_c_loop_int` (interior cells only) instead of `begin_c_loop`,
which also visits the exterior (ghost) cells that mirror the neighbour
partition.

## 3. Move each node once

In `DEFINE_GRID_MOTION`, a mesh node is shared by several faces and, in
parallel, can be seen on more than one partition. The pair
`NODE_POS_NEED_UPDATE(v)` / `NODE_POS_UPDATED(v)` makes sure it moves
once. Skipping it moves shared nodes twice and the mesh folds.

## Things that are safe without special handling

- `DEFINE_CG_MOTION` and `DEFINE_SDOF_PROPERTIES`: Fluent calls them
  consistently on all processes and uses the result; you only need to
  guard the file writing.
- `DEFINE_SOURCE`, `DEFINE_PROPERTY`: called per cell on the node that
  owns it.
- Reading a small input file in every process (the table reader does
  this): fine on a normal installation where all nodes see the working
  directory.

## Testing for parallel bugs without a cluster

Run the same case with `-t2` on your own machine and compare the log
file with the serial run. Compare with a defined tolerance: small floating-point differences are
expected; material discrepancies need investigation. Do this once for every new UDF.
