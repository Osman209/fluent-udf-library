#!/bin/sh
# Runs every check that can run without Fluent:
#   1. gcc syntax check of every UDF against the mock headers (2D and 3D)
#   2. plain-C tests of the shared wave theory and table reader
#   3. Python cross-checks of the C output against numpy
set -e
cd "$(dirname "$0")/.."
echo "--- 1. syntax check of UDFs (mock udf.h) ---"
for nd in 2 3; do
  for f in motion/*.c wavetank/*.c profiles/*.c; do
    gcc -std=c99 -DND_ND=$nd -Wall -Wno-unused-parameter -Wno-unused-variable \
        -Wno-unused-function -Wno-unused-but-set-variable -Wno-misleading-indentation -fsyntax-only \
        -I tests/mock -I common "$f"
  done
  echo "ND_ND=$nd: all files pass"
done
echo "--- 2. C tests ---"
mkdir -p tests/out && cd tests/out
gcc -std=c99 -Wall -O2 -I ../../common ../test_wave_theory.c -o test_wave_theory -lm
gcc -std=c99 -Wall -O2 -I ../../common ../test_table.c -o test_table -lm
gcc -std=c99 -Wall -O2 -I ../../common ../test_profiles.c -o test_profiles -lm
./test_wave_theory
./test_table
./test_profiles
cd ../..
echo "--- 3. Python cross-checks ---"
python3 validation/wave_theory_check.py tests/out
python3 validation/jonswap_check.py tests/out
python3 validation/motion_reference.py
python3 validation/wave_probe_analysis.py --selftest
python3 validation/sdof_reference.py --selftest
python3 validation/profiles_check.py
echo "ALL TESTS PASSED"
