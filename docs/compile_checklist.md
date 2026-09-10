# Compile checklist

Use **compiled** UDFs for everything in this library. The interpreter does
not support `dynamesh_tools.h`, static arrays of this size, or the file
I/O used for logging.

## Before you build

1. A C compiler that Fluent can find. On Windows, Fluent 2023 and later
   ship a built-in Clang; older versions need Visual Studio and Fluent
   started from the VS command prompt or with the compiler in PATH. On
   Linux gcc is found automatically.
2. Put the `.c` file and the headers from `common/` in the working
   directory, or list the headers under "Header Files" in the Compiled
   UDFs panel. The `.c` files include them as `"udf_common.h"` and
   `"wave_theory.h"` with no path.
3. Comments and strings must be plain ASCII. A non-ASCII character
   (Arabic, accented letters, smart quotes copied from a web page) makes
   the Windows build fail with a code page error. This library is
   ASCII-only for that reason.
4. Edit the parameters at the top of the file. Zone IDs
   (`WALL_ZONE_ID`) come from the Boundary Conditions panel, phase index
   (`WATER_PHASE`) from the order in the Phases panel (primary = 0).

## Errors you will see and what they mean

**"The UDF library you are trying to load is not compiled for parallel
use"** - you built in serial and are running in parallel, or the other
way round. Rebuild in the same mode you run in.

**Error 193 / "%1 is not a valid Win32 application"** - the library was
built 32-bit or with a compiler that does not match the Fluent build.
Rebuild with the compiler Fluent expects (64-bit, same version family).

**"ud_io1.h: No such file or directory"** - appears when the build
directory is broken or when you compile from the wrong working
directory. Delete the `libudf` folder, restart Fluent, build again from
the directory that contains the source.

**"udf.h: No such file"** - the compiler is not seeing the Fluent include
path. Usually the build was started outside Fluent, or Fluent was not
started from an environment where the compiler works. Build from inside
the Fluent panel.

**Segmentation fault / "received a fatal signal"** at run time - the most
common causes, in order:
1. A loop over cells or faces running on the host in parallel. Wrap it
   in `#if !RP_HOST ... #endif`.
2. `Lookup_Thread` returned NULL because the zone ID is wrong. Check it
   before use (the logger does).
3. Accessing a phase thread on a single-phase case, or
   `THREAD_SUB_THREAD` with the wrong phase index.
4. Writing to a file from every process at once.
5. A `DEFINE_PROFILE` hooked to the wrong variable (e.g. a VOF profile
   hooked to velocity), which makes Fluent pass a thread the code does
   not expect.

**"Info: 6DOF: can't compute angular acceleration"** printed every
iteration - see `pitfalls_6dof.md`; it is a case setup issue, not a
compile issue.

**Macro not found / wrong number of arguments** - Fluent changed a macro
between versions, or the macro is not documented at all. Check the UDF
manual for your version. One version change is described in
`pitfalls_6dof.md`. Undocumented macros are a separate trap: Ansys does
not publish their behaviour and does not support questions about them on
the forum, so this library prefers documented macros wherever there is a
choice (see the viscous force in `motion/force_moment_logger.c`).

**Comment ends early / a stream of errors from stdio.h** - a `*` followed
by `/` inside a block comment closes it. Writing `u*/kappa` in a comment
does exactly that and the compiler then reads the rest of the file as
code. Write `ustar / kappa` instead.

## After the build

- "Done." in the console with no "error" lines.
- Load, then check that the function names appear in the hook drop-downs.
- Run one time step and read the console: the wave inlet prints k, L and
  the Ursell number; the irregular inlet prints 4 sqrt(m0); the SDOF
  files print the reference CG. If a start-up message is missing, the
  UDF is not hooked.
