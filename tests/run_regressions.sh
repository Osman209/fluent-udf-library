#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
mkdir -p tests/out/regressions
# Syntax matrix: dimensions x precision x serial/host/compute node.
for nd in 2 3; do
 for prec in '' '-DUDF_TEST_SINGLE'; do
  for role in '-DRP_HOST=0 -DRP_NODE=0' '-DRP_HOST=1 -DRP_NODE=0' '-DRP_HOST=0 -DRP_NODE=1'; do
   for f in motion/*.c wavetank/*.c profiles/*.c; do
    gcc -std=c99 -DND_ND=$nd $prec $role -Werror=incompatible-pointer-types -fsyntax-only -I tests/mock -I common "$f"
   done
  done
  for name in actual_cp actual_force actual_spring invalid_tables; do
   gcc -std=c99 -DND_ND=$nd $prec -I tests/mock -I common "tests/test_$name.c" -lm -o "tests/out/regressions/$name"
   (cd tests/out/regressions && "./$name")
  done
 done
done
for role in '-DRP_HOST=0 -DRP_NODE=0 -DEXPECTED_WRITER=1' '-DRP_HOST=1 -DRP_NODE=0 -DEXPECTED_WRITER=0' '-DRP_HOST=0 -DRP_NODE=1 -DI_AM_NODE_ZERO_P=1 -DEXPECTED_WRITER=1' '-DRP_HOST=0 -DRP_NODE=1 -DI_AM_NODE_ZERO_P=0 -DEXPECTED_WRITER=0'; do
 gcc -std=c99 $role -I tests/mock -I common tests/test_writer.c -lm -o tests/out/regressions/writer
 tests/out/regressions/writer
done
printf '%s\n' 'PASS 204 mock syntax configurations; actual-function regressions; 4 writer roles'
