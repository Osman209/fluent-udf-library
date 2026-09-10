# Contributing

Contributions are welcome: a new UDF, a fix, a better explanation, or a
report that something did not work on your version of Fluent.

## Reporting a problem

Open an issue. The templates ask for the Fluent version, the operating
system, whether you were in serial or parallel, and the console output.
That is usually enough to find the cause without going back and forth.

If you would rather not post publicly, email
mohamedosmannn2999@gmail.com.

## Adding a UDF

The house style, so the files stay consistent:

1. **Parameters at the top**, between the two `user parameters` comment
   lines, with units in the comment on each line.
2. **A header comment** that says what the UDF does, the formula it
   implements with a reference where there is one, where to hook it, and
   anything that is easy to get wrong.
3. **ASCII only.** A non-ASCII character in a comment breaks the Windows
   build with a code page error.
4. **Parallel safe.** Loops over cells or faces inside `#if !RP_HOST`,
   sums reduced with `PRF_GRSUM1`, file writing behind `UDF_IS_WRITER`,
   `PRINCIPAL_FACE_P` on faces, `begin_c_loop_int` on cells. See
   `docs/parallel_safety.md`.
5. **Prefer documented macros.** Ansys does not support questions about
   undocumented macros and their behaviour can change without notice. If
   an undocumented macro is the only way, say so in the comment and give
   the user a way to check the result.
6. **Add a check.** If the UDF implements a formula, put the formula in a
   plain-C helper in `common/` where it can be compiled without Fluent,
   and add a case to `tests/` and a numpy cross-check in `validation/`.
   That is what keeps the library trustworthy.

Then run:

    sh tests/run_tests.sh

## Adding a check to an existing file

Also welcome, and often more useful than a new UDF. If you validated one
of these against an experiment or a published case, a script in
`validation/` that reproduces the comparison is a good contribution.
