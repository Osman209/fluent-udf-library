## What this changes

<!-- One or two sentences. -->

## Checks

- [ ] `sh tests/run_tests.sh` passes
- [ ] The file is ASCII only (no accented letters or smart quotes in comments)
- [ ] Any loop over cells or faces is inside `#if !RP_HOST`
- [ ] Any file writing goes through `UDF_IS_WRITER`
- [ ] The user parameters are at the top of the file, in the same style as the others
- [ ] The header comment says what the UDF does and where to hook it

## If you built it in Fluent

Fluent version, operating system, serial or parallel, and what you ran it
on. Not required, but useful for everyone.
